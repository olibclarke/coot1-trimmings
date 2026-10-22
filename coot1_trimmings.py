"""Oli's Coot 1.x customizations.

Copy this file into Coot's startup script directory and restart Coot.
The main user-editable knobs live near the top of the file under
"User-editable settings".
"""

import math
import ast
import base64
import hashlib
import json
import os
import re
import shutil
import struct
import sys
import tempfile
import time
import traceback
import warnings
from bisect import bisect_left
from datetime import datetime, timezone
from functools import lru_cache
from itertools import product
from contextlib import contextmanager, nullcontext

STARTUP_DIAGNOSTIC_MESSAGES = []
_STARTUP_DIAGNOSTIC_DIALOGS_SHOWN = set()
_GTK_REPOSITORY_WARNINGS_SHOWN = set()
_MODEL_EDIT_EPOCHS = {}


def _mark_model_changed(mol_id):
  _MODEL_EDIT_EPOCHS[mol_id] = _MODEL_EDIT_EPOCHS.get(mol_id, 0) + 1


@contextmanager
def _grouped_model_edit(mol_id):
  """Group native undo entries and restore the caller's backup setting."""
  enabled = coot.backup_state(mol_id)
  if enabled not in (0, 1):
    raise ValueError("Cannot determine the model's backup state.")
  if enabled:
    coot.make_backup(mol_id)
  try:
    coot.turn_off_backup(mol_id)
    yield
  finally:
    if enabled:
      coot.turn_on_backup(mol_id)
    else:
      coot.turn_off_backup(mol_id)
    _mark_model_changed(mol_id)


def _record_startup_diagnostic(level, message):
  prefix = f"coot_trimmings startup {level}: "
  full_message = prefix + message
  print(full_message)
  STARTUP_DIAGNOSTIC_MESSAGES.append((level, message))


def _show_startup_diagnostic_once(key, message):
  if key in _STARTUP_DIAGNOSTIC_DIALOGS_SHOWN:
    return None
  _STARTUP_DIAGNOSTIC_DIALOGS_SHOWN.add(key)
  try:
    coot.info_dialog(message)
  except Exception:
    print("coot_trimmings startup notice:", message)
  return None


try:
  import coot
except Exception as error:
  _record_startup_diagnostic(
    "fatal",
    f"Could not import the coot Python module ({error!r}). This script must be run inside Coot.",
  )
  raise

try:
  import coot_utils
except Exception as error:
  fatal_message = (
    "Could not import coot_utils. This Coot build is missing a required Python helper module, "
    f"so coot_trimmings cannot load.\n\nOriginal error: {error!r}"
  )
  _record_startup_diagnostic("fatal", fatal_message)
  _show_startup_diagnostic_once("missing_coot_utils", fatal_message)
  raise

try:
  import gap
except Exception as error:
  gap = None
  _record_startup_diagnostic(
    "warning",
    f"Could not import gap ({error!r}). Gap-fitting helpers will be disabled.",
  )

try:
  from mutate import three_letter_code2single_letter
except Exception as error:
  _record_startup_diagnostic(
    "warning",
    f"Could not import mutate.three_letter_code2single_letter ({error!r}). "
    "Using a built-in residue-code fallback.",
  )

  _THREE_TO_ONE_FALLBACK = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "MSE": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y",
    "VAL": "V", "ASX": "B", "GLX": "Z", "SEC": "U", "PYL": "O",
  }

  def three_letter_code2single_letter(resname):
    return _THREE_TO_ONE_FALLBACK.get(str(resname or "").upper(), "X")

try:
  import gemmi
except Exception:
  gemmi = None

from coot import *
from coot_utils import *


# ============================================================================
# Coot 1.1 compatibility shims
# ============================================================================
#
# Several bundled Python helpers in current Coot 1.x builds are almost usable
# but miss one or two names at runtime. The small wrappers below keep those
# helpers working without having to rewrite the whole upstream module.
#
if not hasattr(coot_utils, "regularize_zone"):
  coot_utils.regularize_zone = coot.regularize_zone
if not hasattr(coot_utils, "get_view_matrix_element"):
  coot_utils.get_view_matrix_element = coot.get_view_matrix_element

_OPTIONAL_COOT_API_WARNINGS_SHOWN = set()
_COOT_RESIDUE_COMPAT_WARNINGS_SHOWN = set()
_raw_residue_info_py = coot.residue_info_py
_raw_residue_name = globals().get("_raw_residue_name", residue_name)
_raw_residue_exists_qm = globals().get("_raw_residue_exists_qm", residue_exists_qm)
_raw_residue_centre_py = globals().get("_raw_residue_centre_py", residue_centre_py)


def _coot_text_arg(value):
  """Normalize Python values for Coot APIs that expect char const *."""
  if value in (None, False):
    return ""
  if isinstance(value, bytes):
    try:
      return value.decode("utf-8")
    except Exception:
      return value.decode("latin-1", errors="ignore")
  return str(value)


def _coot_int_arg(value):
  """Normalize Python values for Coot APIs that expect an integer residue number."""
  if value is None:
    return None
  if isinstance(value, bool):
    return int(value)
  if isinstance(value, int):
    return value
  if isinstance(value, float):
    return int(value)
  if isinstance(value, str):
    stripped_value = value.strip()
    if not stripped_value:
      return None
    try:
      return int(stripped_value)
    except Exception:
      match = re.search(r"[-+]?\d+", stripped_value)
      if match:
        try:
          return int(match.group(0))
        except Exception:
          return None
      return None
  if isinstance(value, (list, tuple)):
    for item in value:
      normalized_item = _coot_int_arg(item)
      if normalized_item is not None:
        return normalized_item
    return None
  for attr_name in ("num", "resno", "seqnum", "value"):
    if hasattr(value, attr_name):
      normalized_item = _coot_int_arg(getattr(value, attr_name))
      if normalized_item is not None:
        return normalized_item
  try:
    return int(value)
  except Exception:
    return None


def _coot_residue_spec_tuple(chain_id, resno, ins_code):
  """Return a normalized residue-spec tuple suitable for Coot C++ bindings."""
  return (_coot_text_arg(chain_id), _coot_int_arg(resno), _coot_text_arg(ins_code))


def _coot_residue_spec_from_spec(spec):
  """Parse common Coot residue-spec payload variants into (chain, resno, ins_code)."""
  if not isinstance(spec, (list, tuple)):
    return None

  candidate_specs = []
  spec_list = list(spec)
  candidate_specs.append(spec_list)

  if len(spec_list) >= 4:
    candidate_specs.append(spec_list[1:4])

  for item in spec_list[:2]:
    if isinstance(item, (list, tuple)):
      nested_spec = list(item)
      candidate_specs.append(nested_spec)
      if len(nested_spec) >= 4:
        candidate_specs.append(nested_spec[1:4])

  seen_candidates = set()
  for candidate in candidate_specs:
    candidate_key = tuple(repr(value) for value in candidate)
    if candidate_key in seen_candidates:
      continue
    seen_candidates.add(candidate_key)
    try:
      chain_id = residue_spec_to_chain_id(candidate)
      resno = residue_spec_to_res_no(candidate)
      ins_code = residue_spec_to_ins_code(candidate)
    except Exception:
      continue
    if chain_id is False or resno is False or ins_code is False:
      continue
    residue_spec = _coot_residue_spec_tuple(chain_id, resno, ins_code)
    if residue_spec[0] and residue_spec[1] is not None:
      return residue_spec
  return None


def residue_info_py(imol, chain_id, resno, ins_code):
  """Compatibility wrapper around Coot's residue_info_py for stricter builds."""
  normalized_chain_id, normalized_resno, normalized_ins_code = _coot_residue_spec_tuple(
    chain_id, resno, ins_code
  )
  if normalized_resno is None:
    warning_key = ("bad_residue_info_resno", repr(resno))
    if warning_key not in _COOT_RESIDUE_COMPAT_WARNINGS_SHOWN:
      _COOT_RESIDUE_COMPAT_WARNINGS_SHOWN.add(warning_key)
      print(
        "coot_trimmings compatibility warning: could not normalize residue number for residue_info_py:",
        repr((chain_id, resno, ins_code)),
      )
    return []
  try:
    return _raw_residue_info_py(imol, normalized_chain_id, normalized_resno, normalized_ins_code) or []
  except TypeError:
    warning_key = ("residue_info_typeerror", repr((normalized_chain_id, normalized_resno, normalized_ins_code)))
    if warning_key not in _COOT_RESIDUE_COMPAT_WARNINGS_SHOWN:
      _COOT_RESIDUE_COMPAT_WARNINGS_SHOWN.add(warning_key)
      print(
        "coot_trimmings compatibility warning: residue_info_py rejected normalized residue spec:",
        repr((normalized_chain_id, normalized_resno, normalized_ins_code)),
      )
      traceback.print_exc()
    return []


def residue_name(imol, chain_id, resno, ins_code):
  """Compatibility wrapper around residue_name for stricter Python/C++ bindings."""
  normalized_chain_id, normalized_resno, normalized_ins_code = _coot_residue_spec_tuple(
    chain_id, resno, ins_code
  )
  if normalized_resno is None:
    return ""
  try:
    return _raw_residue_name(imol, normalized_chain_id, normalized_resno, normalized_ins_code) or ""
  except TypeError:
    warning_key = ("residue_name_typeerror", repr((normalized_chain_id, normalized_resno, normalized_ins_code)))
    if warning_key not in _COOT_RESIDUE_COMPAT_WARNINGS_SHOWN:
      _COOT_RESIDUE_COMPAT_WARNINGS_SHOWN.add(warning_key)
      print(
        "coot_trimmings compatibility warning: residue_name rejected normalized residue spec:",
        repr((normalized_chain_id, normalized_resno, normalized_ins_code)),
      )
      traceback.print_exc()
    return ""


def residue_exists_qm(imol, chain_id, resno, ins_code):
  """Compatibility wrapper around residue_exists_qm for stricter Python/C++ bindings."""
  normalized_chain_id, normalized_resno, normalized_ins_code = _coot_residue_spec_tuple(
    chain_id, resno, ins_code
  )
  if normalized_resno is None:
    return 0
  try:
    return _raw_residue_exists_qm(imol, normalized_chain_id, normalized_resno, normalized_ins_code)
  except TypeError:
    warning_key = ("residue_exists_typeerror", repr((normalized_chain_id, normalized_resno, normalized_ins_code)))
    if warning_key not in _COOT_RESIDUE_COMPAT_WARNINGS_SHOWN:
      _COOT_RESIDUE_COMPAT_WARNINGS_SHOWN.add(warning_key)
      print(
        "coot_trimmings compatibility warning: residue_exists_qm rejected normalized residue spec:",
        repr((normalized_chain_id, normalized_resno, normalized_ins_code)),
      )
      traceback.print_exc()
    return 0


def residue_centre_py(imol, chain_id, resno, ins_code):
  """Compatibility wrapper around residue_centre_py for stricter Python/C++ bindings."""
  normalized_chain_id, normalized_resno, normalized_ins_code = _coot_residue_spec_tuple(
    chain_id, resno, ins_code
  )
  if normalized_resno is None:
    return None
  try:
    return _raw_residue_centre_py(imol, normalized_chain_id, normalized_resno, normalized_ins_code)
  except TypeError:
    warning_key = ("residue_centre_typeerror", repr((normalized_chain_id, normalized_resno, normalized_ins_code)))
    if warning_key not in _COOT_RESIDUE_COMPAT_WARNINGS_SHOWN:
      _COOT_RESIDUE_COMPAT_WARNINGS_SHOWN.add(warning_key)
      print(
        "coot_trimmings compatibility warning: residue_centre_py rejected normalized residue spec:",
        repr((normalized_chain_id, normalized_resno, normalized_ins_code)),
      )
      traceback.print_exc()
    return None


def _call_optional_coot_api(api_name, *args):
  """Call a non-essential Coot API if this build exposes it."""
  api_function = globals().get(api_name, getattr(coot, api_name, None))
  if not callable(api_function):
    if api_name not in _OPTIONAL_COOT_API_WARNINGS_SHOWN:
      _OPTIONAL_COOT_API_WARNINGS_SHOWN.add(api_name)
      _record_startup_diagnostic(
        "warning",
        f"This Coot build does not expose optional API {api_name}(), so that startup tweak was skipped.",
      )
    return None
  return api_function(*args)


def _coot_gui_repository_module(module_name):
  """Use Gtk/Gio/GLib objects exposed by coot_gui only.

  Importing gi.repository directly inside some bundled Coot builds can pull in
  a second GTK/Rsvg stack and break themed icons. Stay on the objects Coot has
  already initialized, and degrade gracefully if they are unavailable.
  """
  existing_module = globals().get(module_name)
  if existing_module is not None:
    return existing_module

  # PyGObject repository modules are normally already loaded by Coot/GTK even
  # when coot_gui does not re-export them as attributes.  Prefer those exact
  # module objects so we stay on Coot's existing GI/GTK stack.
  for module_key in (f"gi.repository.{module_name}", module_name):
    existing_module = sys.modules.get(module_key)
    if existing_module is not None:
      globals()[module_name] = existing_module
      return existing_module

  coot_gui_module = _coot_gui_module()
  if coot_gui_module is not None and hasattr(coot_gui_module, module_name):
    existing_module = getattr(coot_gui_module, module_name)
    globals()[module_name] = existing_module
    return existing_module

  # If GTK is already running under PyGObject, importing another repository
  # namespace (not another GTK version) is safe and is the normal way to obtain
  # GLib.  This fixes builds where coot_gui exposes Gtk/Gio but not GLib.
  if "gi" in sys.modules:
    try:
      repository = __import__(f"gi.repository.{module_name}", fromlist=[module_name])
      globals()[module_name] = repository
      return repository
    except Exception as error:
      _record_startup_diagnostic(
        "warning",
        f"Could not obtain gi.repository.{module_name}: {error}",
      )

  if module_name not in _GTK_REPOSITORY_WARNINGS_SHOWN:
    _GTK_REPOSITORY_WARNINGS_SHOWN.add(module_name)
    _record_startup_diagnostic(
      "warning",
      "This Coot build does not expose coot_gui.{module_name}; "
      "Gtk helper features depending on it will be limited.".format(module_name=module_name),
    )
  return None


def _gap_residue_info_compat(imol, chain_id, resno, ins_code):
  """Match gap.py's older expectation that residue_info() returns a list."""
  return residue_info_py(imol, chain_id, resno, ins_code) or []

if gap is not None:
  gap.residue_info = _gap_residue_info_compat
  if not hasattr(gap, "sys"):
    gap.sys = sys


# ============================================================================
# User-editable settings
# ============================================================================
#
# This section is the main place to tweak script behaviour. If you want to
# change startup defaults, map appearance, lighting behaviour, or a few other
# everyday parameters, start here.
#

# Startup behaviour and display defaults.
STARTUP_SYMMETRY_COLOUR = (255, 35, 0)
STARTUP_USE_LEFT_MOUSE_FOR_ROTATION = 1
STARTUP_SCROLL_BY_WHEEL_MOUSE = 0
STARTUP_AUTO_FIT_BEST_ROTAMER_CLASH = 1
STARTUP_REFINE_MAX_RESIDUES = 100
STARTUP_ROTAMER_SEARCH_MODE = ROTAMERSEARCHLOWRES
STARTUP_MAP_SAMPLING_RATE = 3.0
STARTUP_ADD_TERMINAL_PHI_PSI_TRIALS = 1000
STARTUP_REFINEMENT_WEIGHT = 20.0
STARTUP_SMOOTH_SCROLL = 0
STARTUP_ACTIVE_MAP_DRAG = 0
STARTUP_SHOW_ENVIRONMENT_DISTANCES = 0
STARTUP_SHOW_ENVIRONMENT_BUMPS = 0
STARTUP_SHOW_ENVIRONMENT_H_BONDS = 1
STARTUP_ENVIRONMENT_DISTANCE_LIMITS = (2.1, 3.2)
STARTUP_MAP_RADIUS = 20
STARTUP_MAP_RADIUS_EM = 20
STARTUP_DEFAULT_BOND_THICKNESS = 3
STARTUP_USE_VARIABLE_BOND_THICKNESS = 1
STARTUP_DEFAULT_NEW_ATOM_B = 50.0
STARTUP_ADD_TERMINAL_POST_REFINE = 1
STARTUP_TERMINAL_RIGID_BODY_REFINE = 0
STARTUP_MUTATE_AUTO_FIT_POST_REFINE = 1
STARTUP_SYMMETRY_SIZE = 30
STARTUP_NOMENCLATURE_ERRORS_MODE = "ignore"
STARTUP_REORIENTING_NEXT_RESIDUE_MODE = 1
STARTUP_CURMUDGEON_MODE = 1
STARTUP_ROTATION_CENTRE_CROSSHAIRS_SCALE = 0.15
STARTUP_ROTATION_CENTRE_CROSSHAIRS_COLOUR = (1.0, 1.0, 1.0, 1.0)

# User-facing helper defaults.
DEFAULT_NEW_HELIX_CHAIN_ID = "A"
HIGH_CONTRAST_BOND_THICKNESS = 2
PROPORTIONAL_EDITING_RADIUS = 0.0
PROPORTIONAL_EDITING_RADIUS_STATE = {"applied": 0.0, "pending": None, "source": None}
COMMON_MONOMER_FAVORITES_FILENAME = "coot_trimmings_favorites.json"
DISPLAY_PREFERENCES_FILENAME = "coot_trimmings_preferences.py"
CUSTOM_TOOLBAR_FILENAME = "coot_trimmings_toolbar.json"
TRACKPAD_PREFERENCES_FILENAME = "coot_trimmings_trackpad.json"
COMMON_MONOMER_FAVORITE_CIF_PREFIX = "coot_trimmings_favorite_"
WATER_REFINE_MAX_DISPLACEMENT = 1.0
RESIDUE_ANNOTATION_DEFAULT_AUTHOR = (
  os.environ.get("USER")
  or os.environ.get("USERNAME")
  or ""
)
RESIDUE_ANNOTATION_SCHEMA_VERSION = "4"
RESIDUE_ANNOTATION_NEARBY_RADIUS = 6.0
RESIDUE_ANNOTATION_POLL_INTERVAL_MS = 300

# Map restyling defaults for the EM helper.
EM_REFINED_MAP_COLOUR = (26 / 255.0, 95 / 255.0, 180 / 255.0)  # Ambient-map dark blue: #1a5fb4.
EM_MAP_AMBIENT = (1.0, 1.0, 1.0, 1.0)
EM_MAP_DIFFUSE = (0.0, 0.0, 0.0, 1.0)
EM_MAP_SPECULAR = (0.0, 64.0)
EM_MAP_FRESNEL = (0, 0.0, 0.3, 3.0)
EM_REFINED_MAP_CONTOUR_SIGMA = 2.3
EM_TARGET_PIXEL_SIZE = 0.5
EM_RESAMPLE_CONFIRM_MAX_GRID_DIMENSION = 512
EM_MASKED_MAP_RADIUS = 1.7
EM_MASKED_MAP_COLOUR = (1.0, 0.0, 1.0)
MAP_BRIGHTEN_SCALE_FACTOR = 1.05
MAP_DARKEN_SCALE_FACTOR = 1.0 / MAP_BRIGHTEN_SCALE_FACTOR
MAP_BRIGHTNESS_MIN_COMPONENT = 0.05
MAP_BRIGHTNESS_MAX_COMPONENT = 1.0

# EM-ringer-like side-chain density helper defaults.
EMRINGER_HELPER_STEP_DEGREES = 15.0
EMRINGER_HELPER_OCCUPIED_DISTANCE = 1.2
EMRINGER_HELPER_CA_CB_BOND_LENGTH = 1.53
EMRINGER_HELPER_CB_CG_BOND_LENGTH = 1.52
EMRINGER_HELPER_FINE_STEP_DEGREES = 5.0
EMRINGER_HELPER_MIN_BACKBONE_SUPPORT_FRACTION = 0.75
EMRINGER_HELPER_WEAK_LIGAND_OUTSIDE_FRACTION = 0.30
EMRINGER_HELPER_WEAK_STAGE_INSET_BUFFER = 0.30
EMRINGER_HELPER_MISFIT_MIN_PEAK_TO_CONTOUR_RATIO = 1.25
EMRINGER_HELPER_MISFIT_MIN_GAIN_TO_CONTOUR_RATIO = 0.50
EMRINGER_HELPER_MISFIT_MIN_PROMINENCE_TO_CONTOUR_RATIO = 0.25
EMRINGER_HELPER_MISFIT_SUPPORT_STEP_DEGREES = 5.0

# Thresholded Q-residue contour-optimization prototype defaults.
THRESHOLDED_Q_GAUSSIAN_SIGMA = 0.60
THRESHOLDED_Q_IDEAL_ENCLOSED_MASS_PERCENT = 70.0
THRESHOLDED_Q_GRID_SPACING = 0.60
THRESHOLDED_Q_SAMPLE_RADIUS_SIGMA = 4.0
THRESHOLDED_Q_SAMPLE_SHELL_PADDING = 0.75
THRESHOLDED_Q_CONTEXT_ATOM_INFLUENCE_SIGMA = 4.0
THRESHOLDED_Q_LOCAL_RADIUS = 6.0
THRESHOLDED_Q_CONTEXT_RADIUS_PADDING = 4.0
THRESHOLDED_Q_LOCAL_MAX_RESIDUES = 12
THRESHOLDED_Q_SWEEP_MIN_SIGMA = 0.5
THRESHOLDED_Q_SWEEP_MAX_SIGMA = 12.0
THRESHOLDED_Q_SWEEP_STEPS = 50

_THRESHOLDED_Q_HYDROGEN_ELEMENTS = {"H", "D"}
_THRESHOLDED_Q_BACKBONE_ATOM_NAMES = {"N", "CA", "C", "O", "P", "OP1", "OP2", "O1P", "O2P", "O3P"}

# Model lighting presets used by the high-contrast toggle.
MODEL_NORMAL_AMBIENT = (0.35, 0.35, 0.35, 1.0)
MODEL_NORMAL_DIFFUSE = (0.75, 0.75, 0.75, 1.0)
MODEL_NORMAL_SPECULAR = (0.20, 64.0)
MODEL_AMBIENT_ONLY_AMBIENT = (1.0, 1.0, 1.0, 1.0)
MODEL_AMBIENT_ONLY_DIFFUSE = (0.0, 0.0, 0.0, 1.0)
MODEL_AMBIENT_ONLY_SPECULAR = (0.0, 64.0)


# ============================================================================
# Runtime state
# ============================================================================

REGISTERED_KEYBINDING_CALLBACKS = []
REGISTERED_KEYBINDINGS = {}
REGISTERED_KEYBINDING_ACTIONS = {}
_CUSTOM_KEY_BINDINGS_WINDOW = None
_CUSTOM_TOOLBAR_STATE = {}
_TRACKPAD_EDIT_STATE = None
_TRACKPAD_STATE = {"enabled": False, "bindings": [], "native_zoom": [], "native_scroll": [], "claimed_scroll": {}, "pending": {}, "source": None,
                   "pinch": {}, "scroll_active": {}, "window": None, "reverse": True, "startup": False,
                   "save_source": None, "startup_source": None, "startup_attempts": 0,
                   "sensitivity": {"rotate": 2.0, "translate": 2.0, "zoom": 2.0, "contour": 2.0, "clipping": 2.0}}
MAP_SURFACE_DISPLAY_STATE = {}
MAP_SURFACE_OPACITY_STATE = {}
MAP_GLOBAL_VIEW_SETTINGS = {}
MAP_LOCAL_APPEARANCE_STATE = {}
SMART_COPY_TEMPLATE_IMOL = None
SMART_COPY_SOURCE_CENTRE = None
SMART_COPY_RESIDUE_NAME = None
MODEL_AMBIENT_LIGHTING_ENABLED = False
REFINEMENT_RAMACHANDRAN_OVERLAY_ENABLED = globals().get("REFINEMENT_RAMACHANDRAN_OVERLAY_ENABLED", True)
MODEL_PRE_HIGH_CONTRAST_GL_LIGHTING_STATE = None
MODEL_PRE_HIGH_CONTRAST_BOND_THICKNESS = None
MODEL_HIGH_CONTRAST_MOLECULES = set()
NAVIGATION_LAST_RESIDUE = None
COMMON_MONOMER_FAVORITES_MENU = None
COMMON_MONOMER_LAST_BROWSED_DIRECTORY = None
RESIDUE_ANNOTATION_LAST_AUTHOR = RESIDUE_ANNOTATION_DEFAULT_AUTHOR
RESIDUE_ANNOTATIONS_BY_MOLECULE = {}
EMRINGER_HELPER_BACKBONE_ATOM_NAMES = {"N", "CA", "C", "O", "OXT", "CB"}
RESIDUE_ANNOTATION_CATEGORY_PREFIX = "_cootnote_residue_note."
RESIDUE_ANNOTATION_META_VERSION_TAG = "_cootnote_annotation_meta.version"
RESIDUE_ANNOTATION_META_PROGRAM_TAG = "_cootnote_annotation_meta.program"
RESIDUE_ANNOTATION_RANGE_FIELDS = (
  "target_kind", "end_label_comp_id", "end_label_asym_id", "end_label_seq_id",
  "end_auth_asym_id", "end_auth_seq_id", "end_pdbx_PDB_ins_code",
)
RESIDUE_ANNOTATION_REGION_HIGHLIGHTS = {}
RESIDUE_ANNOTATION_RANGE_PICK_STATE = {}


@lru_cache(maxsize=1)
def _coot_trimmings_file_path():
  candidates = []

  runtime_file = globals().get("__file__", "")
  module_file = getattr(sys.modules.get(__name__), "__file__", "")

  for candidate in (
    runtime_file,
    module_file,
    os.path.expanduser("~/.config/Coot/coot_trimmings.py"),
    os.path.expanduser("~/Library/Application Support/Coot/coot_trimmings.py"),
    os.path.join(os.getcwd(), "coot_trimmings.py"),
  ):
    if isinstance(candidate, str) and candidate:
      normalized = os.path.abspath(candidate)
      if normalized not in candidates:
        candidates.append(normalized)

  for candidate in candidates:
    if os.path.isfile(candidate):
      return candidate

  if candidates:
    return candidates[0]
  return os.path.abspath("coot_trimmings.py")


@lru_cache(maxsize=1)
def _coot_trimmings_dir():
  return os.path.dirname(_coot_trimmings_file_path())



def fit_gap(imol, chain_id, start_resno, stop_resno, sequence="", use_rama_restraints=1):
  """Coot 1.x-safe wrapper around the bundled gap fitter.

  The stock Python helper in this build has a few runtime bugs and leaves
  temporary molecules behind. This wrapper keeps the original "fit both
  directions and keep the better answer" logic, while making cleanup and error
  handling predictable.
  """
  if gap is None:
    add_status_bar_text("Gap fitting is unavailable in this Coot Python build")
    return None
  imol_map = coot.imol_refinement_map()
  if imol_map == -1:
    coot.info_dialog("Need to set a map to fit a loop")
    return None

  backup_mode = coot.backup_state(imol)
  coot.make_backup(imol)
  coot.turn_off_backup(imol)

  rama_status = coot.refine_ramachandran_angles_state()
  coot.set_refine_ramachandran_angles(use_rama_restraints)
  temp_imols = set()

  try:
    if stop_resno < start_resno:
      res_limits = [stop_resno - 1, start_resno + 1]
    else:
      res_limits = [start_resno - 1, stop_resno + 1]

    if all([residue_exists_qm(imol, chain_id, resno, "") for resno in res_limits]):
      existing_models = set(model_molecule_list())
      imol_backwards = coot.copy_molecule(imol)
      if imol_backwards in existing_models or not valid_model_molecule_qm(imol_backwards):
        raise RuntimeError("Could not create an independent backwards fitting model.")
      temp_imols.add(imol_backwards)
      loop_len = abs(start_resno - stop_resno) + 1
      imol_both = coot.copy_molecule(imol) if loop_len >= 6 else None
      if imol_both is not None:
        if imol_both in existing_models or imol_both == imol_backwards or not valid_model_molecule_qm(imol_both):
          raise RuntimeError("Could not create an independent combined fitting model.")
        temp_imols.add(imol_both)

      atom_selection = "//" + chain_id + "/" + str(min(start_resno, stop_resno) - 1) + \
                       "-" + str(max(start_resno, stop_resno) + 1)
      imol_fragment_backup = coot.new_molecule_by_atom_selection(imol, atom_selection)
      if imol_fragment_backup in existing_models or imol_fragment_backup in temp_imols or not valid_model_molecule_qm(imol_fragment_backup):
        raise RuntimeError("Could not create an independent gap backup fragment.")
      temp_imols.add(imol_fragment_backup)
      coot.set_mol_displayed(imol_fragment_backup, 0)

      gap.fit_gap_generic(imol, chain_id, start_resno, stop_resno, sequence)
      gap.fit_gap_generic(imol_backwards, chain_id, stop_resno, start_resno, sequence)

      result_a = gap.low_density_average(imol_map, imol, chain_id, start_resno, stop_resno)
      result_b = gap.low_density_average(imol_map, imol_backwards, chain_id, start_resno, stop_resno)
      loop_list = [[imol, result_a], [imol_backwards, result_b]]

      if loop_len >= 6:
        start_resno1 = min([start_resno, stop_resno])
        stop_resno2 = max([start_resno, stop_resno])
        half_loop_len = loop_len // 2
        stop_resno1 = start_resno1 + half_loop_len - 1
        start_resno2 = stop_resno1 + 1
        sequence1 = sequence[0:half_loop_len]
        sequence2 = sequence[half_loop_len:len(sequence)]
        gap.fit_gap_generic(imol_both, chain_id, start_resno1, stop_resno1, sequence1)
        gap.fit_gap_generic(imol_both, chain_id, stop_resno2, start_resno2, sequence2)
        immediate_refinement_mode = coot.refinement_immediate_replacement_state()
        try:
          coot.set_refinement_immediate_replacement(1)
          coot.refine_zone(imol_both, chain_id, stop_resno1 - loop_len//3, start_resno2 + loop_len//3, "")
          coot.accept_regularizement()
        finally:
          coot.set_refinement_immediate_replacement(immediate_refinement_mode)
        result_c = gap.low_density_average(imol_map, imol_both, chain_id, start_resno, stop_resno)
        loop_list.append([imol_both, result_c])

      candidates = [(candidate, float(score)) for candidate, score in loop_list
                    if valid_model_molecule_qm(candidate) and math.isfinite(float(score))]
      if not candidates:
        raise RuntimeError("Gap fitting produced no valid scored candidate.")
      winner, _score = max(candidates, key=lambda item: item[1])
      if winner != imol and not coot.replace_fragment(imol, winner, atom_selection):
        raise RuntimeError("Coot rejected the best gap-fit fragment.")
    else:
      gap.fit_gap_generic(imol, chain_id, start_resno, stop_resno, sequence)
  finally:
    for temp_imol in list(temp_imols):
      if valid_model_molecule_qm(temp_imol):
        coot.close_molecule(temp_imol)
    coot.set_refine_ramachandran_angles(rama_status)
    if backup_mode == 1:
      coot.turn_on_backup(imol)
    _mark_model_changed(imol)

if gap is not None:
  gap.fit_gap = fit_gap


def _residue_spec_to_cid_compat(residue_spec):
  """Convert a simple residue spec into the CID style used by newer helpers."""
  if not isinstance(residue_spec, (list, tuple)) or len(residue_spec) < 3:
    return None
  chain_id = residue_spec[0]
  res_no = residue_spec[1]
  ins_code = residue_spec[2] or ""
  if ins_code:
    return "//{chain_id}/{res_no}.{ins_code}".format(
      chain_id=chain_id,
      res_no=res_no,
      ins_code=ins_code,
    )
  return "//{chain_id}/{res_no}".format(chain_id=chain_id, res_no=res_no)


USER_DEFINED_COLOUR_TABLE_SIZE = 60
USER_DEFINED_COLOUR_TABLE_BASE = 60
USER_DEFINED_COLOUR_TABLE_READY = False
USER_DEFINED_DIRECT_COLOUR_MIN_REVISION_COUNT = 19935
USER_DEFINED_DIRECT_COLOUR_MIN_COMMIT_PREFIX = "bcc3af8"


def _make_default_user_defined_colours():
  explicit_colours = {
    0:  [0.72, 0.72, 0.72],
    1:  [0.92, 0.92, 0.92],
    2:  [0.18, 0.38, 0.95],
    3:  [0.20, 0.62, 0.98],
    4:  [0.10, 0.22, 0.98],
    5:  [0.42, 0.78, 1.00],
    10: [0.00, 0.76, 0.68],
    15: [0.10, 0.78, 0.18],
    22: [0.98, 0.88, 0.18],
    27: [0.98, 0.60, 0.12],
    28: [1.00, 0.48, 0.02],
    30: [1.00, 0.28, 0.08],
    31: [0.92, 0.12, 0.12],
    34: [0.90, 0.20, 0.85],
    39: [0.55, 0.28, 0.90],
    # Density-fit palette: ordered from poor (red) to good (blue), avoiding
    # magenta so the meaning stays visually obvious.
    48: [0.78, 0.00, 0.12],
    49: [0.90, 0.24, 0.08],
    50: [0.98, 0.50, 0.06],
    51: [1.00, 0.74, 0.10],
    52: [0.96, 0.90, 0.20],
    53: [0.72, 0.88, 0.18],
    54: [0.34, 0.86, 0.22],
    55: [0.10, 0.78, 0.42],
    56: [0.06, 0.78, 0.68],
    57: [0.08, 0.70, 0.90],
    58: [0.12, 0.52, 0.98],
    59: [0.18, 0.36, 0.98],
  }

  def interpolated_colour(colour_index):
    if colour_index in explicit_colours:
      return list(explicit_colours[colour_index])

    lower_indices = [idx for idx in explicit_colours if idx < colour_index]
    upper_indices = [idx for idx in explicit_colours if idx > colour_index]

    if not lower_indices and not upper_indices:
      return [0.72, 0.72, 0.72]
    if not lower_indices:
      return list(explicit_colours[min(upper_indices)])
    if not upper_indices:
      return list(explicit_colours[max(lower_indices)])

    lower_index = max(lower_indices)
    upper_index = min(upper_indices)
    lower_colour = explicit_colours[lower_index]
    upper_colour = explicit_colours[upper_index]
    blend = float(colour_index - lower_index) / float(upper_index - lower_index)
    return [
      lower_colour[i] + blend * (upper_colour[i] - lower_colour[i])
      for i in range(3)
    ]

  base_colours = []
  for colour_index in range(USER_DEFINED_COLOUR_TABLE_SIZE):
    base_colours.append(interpolated_colour(colour_index))

  # CA user-colour rendering indexes this list by bond colour number, whereas
  # all-atom rendering looks up the explicit key. Pad to align positions with
  # our 60..119 keys, repeating the neutral key so standard 0..59 colours are
  # never overridden. Coot retains these duplicate entries in order.
  colours = [
    (USER_DEFINED_COLOUR_TABLE_BASE, list(base_colours[0]))
    for _ in range(USER_DEFINED_COLOUR_TABLE_BASE)
  ]
  for i in range(USER_DEFINED_COLOUR_TABLE_SIZE):
    colours.append((USER_DEFINED_COLOUR_TABLE_BASE + i, base_colours[i]))
  return colours


def _legacy_user_colour_index_to_coot_index(colour_index):
  if not isinstance(colour_index, int):
    return colour_index
  if colour_index < 0:
    return USER_DEFINED_COLOUR_TABLE_BASE
  if colour_index >= USER_DEFINED_COLOUR_TABLE_SIZE:
    return colour_index
  return USER_DEFINED_COLOUR_TABLE_BASE + colour_index


DENSITY_FIT_COLOUR_INDICES = [48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59]
DENSITY_FIT_COLOUR_BAND_NAMES = ["red", "orange", "yellow", "green", "cyan", "blue"]


def _density_fit_score_to_colour_index(score):
  """Map CC score 0..1 to an ordered bad->good colour ramp."""
  if score < 0.0:
    score = 0.0
  if score > 1.0:
    score = 1.0
  palette_size = len(DENSITY_FIT_COLOUR_INDICES)
  palette_index = int(score * palette_size)
  if palette_index >= palette_size:
    palette_index = palette_size - 1
  return DENSITY_FIT_COLOUR_INDICES[palette_index]


def _colour_band_legend_text(metric_label, min_value, max_value, reverse=False):
  band_count = len(DENSITY_FIT_COLOUR_BAND_NAMES)
  span = float(max_value) - float(min_value)
  lines = []
  for band_index, colour_name in enumerate(DENSITY_FIT_COLOUR_BAND_NAMES):
    if reverse:
      low = min_value + span * float(band_count - band_index - 1) / float(band_count)
      high = min_value + span * float(band_count - band_index) / float(band_count)
    else:
      low = min_value + span * float(band_index) / float(band_count)
      high = min_value + span * float(band_index + 1) / float(band_count)
    if band_index == band_count - 1 and not reverse:
      high_text = "{:.2f}".format(max_value)
    elif band_index == 0 and reverse:
      high_text = "{:.2f}".format(max_value)
    else:
      high_text = "{:.2f}".format(high)
    lines.append(
      "{colour}: {metric} {low:.2f}-{high}".format(
        colour=colour_name,
        metric=metric_label,
        low=low,
        high=high_text,
      )
    )
  return "\n".join(lines)


def ensure_user_defined_colour_table():
  global USER_DEFINED_COLOUR_TABLE_READY
  # Always re-apply our legacy palette here. Other helpers (for example
  # Alphafold colouring utilities shipped with Coot) can mutate the global
  # user-defined colour table, so a one-time initialization is not reliable.
  colours = _make_default_user_defined_colours()
  multimap = getattr(coot, "_trimmings_multimap_refinement", None)
  if multimap is not None:
    colours = _make_multimap_user_defined_colours(multimap)
  coot.set_user_defined_colours_py(colours)
  USER_DEFINED_COLOUR_TABLE_READY = True


def _coot_git_commit_string():
  git_commit_function = globals().get("git_commit", getattr(coot, "git_commit", None))
  if not callable(git_commit_function):
    return ""
  try:
    return str(git_commit_function() or "").strip()
  except Exception:
    return ""


def _coot_revision_count():
  version_extra_function = globals().get(
    "coot_version_extra_info",
    getattr(coot, "coot_version_extra_info", None),
  )
  if not callable(version_extra_function):
    return None
  try:
    version_extra = str(version_extra_function() or "")
  except Exception:
    return None
  match = re.search(r"revision-count\s+(\d+)", version_extra)
  if not match:
    return None
  try:
    return int(match.group(1))
  except Exception:
    return None


def _supports_direct_user_defined_colouring():
  revision_count = _coot_revision_count()
  if revision_count is not None and revision_count >= USER_DEFINED_DIRECT_COLOUR_MIN_REVISION_COUNT:
    return True
  git_commit_string = _coot_git_commit_string().lower()
  if git_commit_string.startswith(USER_DEFINED_DIRECT_COLOUR_MIN_COMMIT_PREFIX):
    return True
  return False


if "set_user_defined_atom_colour_by_residue_py" not in globals():
  def set_user_defined_atom_colour_by_residue_py(mol_id, residue_specs_colour_index_tuple_list_py):
    ensure_user_defined_colour_table()
    selection_colour_list = []
    for item in residue_specs_colour_index_tuple_list_py or []:
      if not isinstance(item, (list, tuple)) or len(item) < 2:
        continue
      residue_spec = item[0]
      colour_index = item[1]
      cid = _residue_spec_to_cid_compat(residue_spec)
      if cid is None:
        continue
      selection_colour_list.append((cid, _legacy_user_colour_index_to_coot_index(colour_index)))
    return coot.set_user_defined_atom_colour_by_selection_py(mol_id, selection_colour_list)


def _set_user_defined_atom_colour_by_residue_rows_py(mol_id, residue_specs_colour_index_tuple_list_py, palette_ready=False):
  if not palette_ready:
    ensure_user_defined_colour_table()
  selection_colour_list = []
  for item in residue_specs_colour_index_tuple_list_py or []:
    if not isinstance(item, (list, tuple)) or len(item) < 2:
      continue
    residue_spec = item[0]
    colour_index = item[1]
    cid = residue_spec if isinstance(residue_spec, str) else _residue_spec_to_cid_compat(residue_spec)
    if cid is None:
      continue
    selection_colour_list.append((cid, _legacy_user_colour_index_to_coot_index(colour_index)))
  if not selection_colour_list:
    return 0
  try:
    coot.set_user_defined_atom_colour_by_selection_py(mol_id, selection_colour_list)
  except Exception:
    traceback.print_exc()
    return 0
  return len(selection_colour_list)


if "clear_user_defined_atom_colours" not in globals():
  clear_user_defined_atom_colours = coot.clear_user_defined_atom_colours


if "graphics_to_user_defined_atom_colours_representation" not in globals():
  graphics_to_user_defined_atom_colours_representation = (
    coot.graphics_to_user_defined_atom_colours_representation
  )


if "graphics_to_user_defined_atom_colours_all_atoms_representation" not in globals():
  graphics_to_user_defined_atom_colours_all_atoms_representation = (
    coot.graphics_to_user_defined_atom_colours_all_atoms_representation
  )


def _native_representation_function(name):
  # Keep originals on the module: startup scripts may be re-executed in place.
  originals = getattr(coot, "_trimmings_native_representations", None)
  if originals is None:
    originals = {}
    coot._trimmings_native_representations = originals
  if name not in originals:
    originals[name] = getattr(coot, name)
  return originals[name]


_coot_graphics_to_bonds_representation = _native_representation_function("graphics_to_bonds_representation")
_coot_graphics_to_rainbow_representation = _native_representation_function("graphics_to_rainbow_representation")
_coot_graphics_to_b_factor_representation = _native_representation_function("graphics_to_b_factor_representation")
_coot_graphics_to_ca_plus_ligands_representation = _native_representation_function("graphics_to_ca_plus_ligands_representation")
_coot_graphics_to_ca_plus_ligands_and_sidechains_representation = _native_representation_function("graphics_to_ca_plus_ligands_and_sidechains_representation")
_coot_graphics_to_ca_plus_ligands_sec_struct_representation = _native_representation_function("graphics_to_ca_plus_ligands_sec_struct_representation")
LAST_STANDARD_REPRESENTATION_BY_MOL = {}


def _clear_user_defined_colours_for_standard_representation(mol_id):
  if mol_id is None or mol_id < 0:
    return
  if _supports_direct_user_defined_colouring():
    return
  clear_user_defined_atom_colours(mol_id)


def _remember_standard_representation(mol_id, representation_function):
  if mol_id is None or mol_id < 0:
    return
  LAST_STANDARD_REPRESENTATION_BY_MOL[mol_id] = representation_function
  sequence_function = globals().get("_representation_sequence")
  if callable(sequence_function):
    sequence = sequence_function()
    if representation_function in sequence:
      cycle_rep_flag[mol_id] = sequence.index(representation_function)


def _last_standard_representation_function(mol_id):
  return LAST_STANDARD_REPRESENTATION_BY_MOL.get(
    mol_id,
    graphics_to_bonds_representation,
  )


def _restore_last_standard_representation(mol_id, representation_function=None):
  if representation_function is None:
    representation_function = _last_standard_representation_function(mol_id)
  try:
    return representation_function(mol_id)
  except Exception:
    return _coot_graphics_to_bonds_representation(mol_id)


def graphics_to_bonds_representation(mol_id):
  _clear_user_defined_colours_for_standard_representation(mol_id)
  if _supports_direct_user_defined_colouring():
    get_bond_type = getattr(coot, "get_graphics_molecule_bond_type", None)
    if callable(get_bond_type) and get_bond_type(mol_id) in (12, 13):
      # Leave the custom bond mode before rebuilding normal element colours.
      _coot_graphics_to_ca_plus_ligands_representation(mol_id)
  result = _coot_graphics_to_bonds_representation(mol_id)
  _remember_standard_representation(mol_id, graphics_to_bonds_representation)
  return result


def graphics_to_rainbow_representation(mol_id):
  _clear_user_defined_colours_for_standard_representation(mol_id)
  result = _coot_graphics_to_rainbow_representation(mol_id)
  _remember_standard_representation(mol_id, graphics_to_rainbow_representation)
  return result


def graphics_to_b_factor_representation(mol_id):
  _clear_user_defined_colours_for_standard_representation(mol_id)
  result = _coot_graphics_to_b_factor_representation(mol_id)
  _remember_standard_representation(mol_id, graphics_to_b_factor_representation)
  return result


def graphics_to_ca_plus_ligands_representation(mol_id):
  _clear_user_defined_colours_for_standard_representation(mol_id)
  result = _coot_graphics_to_ca_plus_ligands_representation(mol_id)
  _remember_standard_representation(mol_id, graphics_to_ca_plus_ligands_representation)
  return result


def graphics_to_ca_plus_ligands_and_sidechains_representation(mol_id):
  _clear_user_defined_colours_for_standard_representation(mol_id)
  result = _coot_graphics_to_ca_plus_ligands_and_sidechains_representation(mol_id)
  _remember_standard_representation(mol_id, graphics_to_ca_plus_ligands_and_sidechains_representation)
  return result


def graphics_to_ca_plus_ligands_sec_struct_representation(mol_id):
  _clear_user_defined_colours_for_standard_representation(mol_id)
  result = _coot_graphics_to_ca_plus_ligands_sec_struct_representation(mol_id)
  _remember_standard_representation(mol_id, graphics_to_ca_plus_ligands_sec_struct_representation)
  return result


coot.graphics_to_bonds_representation = graphics_to_bonds_representation
coot.graphics_to_rainbow_representation = graphics_to_rainbow_representation
coot.graphics_to_b_factor_representation = graphics_to_b_factor_representation
coot.graphics_to_ca_plus_ligands_representation = graphics_to_ca_plus_ligands_representation
coot.graphics_to_ca_plus_ligands_and_sidechains_representation = (
  graphics_to_ca_plus_ligands_and_sidechains_representation
)
coot.graphics_to_ca_plus_ligands_sec_struct_representation = (
  graphics_to_ca_plus_ligands_sec_struct_representation
)
CUSTOM_COLOUR_ADDITIONAL_REPRESENTATIONS = {}
CUSTOM_COLOUR_ADDITIONAL_REPRESENTATION_TYPE = 1
CUSTOM_COLOUR_USER_DEFINED_BONDS_BOX_TYPE = 12
CUSTOM_COLOUR_ADDITIONAL_BOND_WIDTH = 0.14

POLYMER_RESIDUE_NAMES = {
  # Standard RNA/DNA.
  "A", "C", "G", "U", "T", "DA", "DC", "DG", "DT", "DU",
  # Common modified nucleotides seen in experimental models.
  "PSU", "H2U", "5MU", "OMU", "4SU",
  "1MA", "2MA", "6MA", "MIA",
  "5MC", "OMC",
  "1MG", "2MG", "7MG", "M2G", "OMG", "YG",
  # Standard amino acids.
  "ALA", "UNK", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN",
  "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "MSE", "PHE",
  "PRO", "SER", "THR", "TRP", "TYR", "VAL",
  # Common modified amino acids / PTMs that should still count as polymer.
  "P1L", "SEP", "TPO", "PTR", "CSO", "CME", "MLY", "MLZ",
  "HYP", "KCX", "ALY", "MLY", "FME", "PYL", "SEC"
}

# Residue-specific atom names that genuinely occupy the chi1/gamma position.
# These are used to decide whether a strong virtual-CG peak is already
# accounted for by placed atoms, without letting distal atoms mask real chi1
# problems.
EMRINGER_HELPER_CHI1_BLOCKER_ATOM_NAMES = {
  "ALA": set(),
  "ARG": {"CG"},
  "ASN": {"CG"},
  "ASP": {"CG"},
  "CME": {"SG"},
  "CSO": {"SG"},
  "CYS": {"SG"},
  "FME": {"CG"},
  "GLN": {"CG"},
  "GLU": {"CG"},
  "HIS": {"CG"},
  "HYP": {"CG"},
  "ILE": {"CG1", "CG2"},
  "KCX": {"CG"},
  "LEU": {"CG"},
  "LYS": {"CG"},
  "M3L": {"CG"},
  "MET": {"CG"},
  "MSE": {"CG"},
  "MLY": {"CG"},
  "MLZ": {"CG"},
  "P1L": {"SG"},
  "PHE": {"CG"},
  "PRO": {"CG"},
  "PTR": {"CG"},
  "SEC": {"SE"},
  "SEP": {"OG"},
  "SER": {"OG"},
  "THR": {"OG1", "CG2"},
  "TPO": {"OG1", "CG2"},
  "TRP": {"CG"},
  "TYR": {"CG"},
  "VAL": {"CG1", "CG2"},
}

# Extra torsion stages for longer flexible sidechains. Each entry defines a
# later sidechain sweep around a downstream bond plus the atom(s) that should
# already occupy that stage if the model is built correctly.
EMRINGER_HELPER_LATE_STAGE_DEFINITIONS = {
  "ARG": [
    {"label": "chi2", "reference": "CA", "axis_start": "CB", "axis_end": "CG", "blockers": {"CD"}, "bond_length": 1.52},
    {"label": "chi3", "reference": "CB", "axis_start": "CG", "axis_end": "CD", "blockers": {"NE"}, "bond_length": 1.46},
    {"label": "chi4", "reference": "CG", "axis_start": "CD", "axis_end": "NE", "blockers": {"CZ"}, "bond_length": 1.33},
  ],
  "GLN": [
    {"label": "chi2", "reference": "CA", "axis_start": "CB", "axis_end": "CG", "blockers": {"CD"}, "bond_length": 1.52},
  ],
  "GLU": [
    {"label": "chi2", "reference": "CA", "axis_start": "CB", "axis_end": "CG", "blockers": {"CD"}, "bond_length": 1.52},
  ],
  "ILE": [
    {"label": "chi2", "reference": "CA", "axis_start": "CB", "axis_end": "CG1", "blockers": {"CD1"}, "bond_length": 1.52},
  ],
  "LEU": [
    {"label": "chi2", "reference": "CA", "axis_start": "CB", "axis_end": "CG", "blockers": {"CD1", "CD2"}, "bond_length": 1.52},
  ],
  "LYS": [
    {"label": "chi2", "reference": "CA", "axis_start": "CB", "axis_end": "CG", "blockers": {"CD"}, "bond_length": 1.52},
    {"label": "chi3", "reference": "CB", "axis_start": "CG", "axis_end": "CD", "blockers": {"CE"}, "bond_length": 1.52},
    {"label": "chi4", "reference": "CG", "axis_start": "CD", "axis_end": "CE", "blockers": {"NZ"}, "bond_length": 1.48},
  ],
  "MET": [
    {"label": "chi2", "reference": "CA", "axis_start": "CB", "axis_end": "CG", "blockers": {"SD"}, "bond_length": 1.81},
    {"label": "chi3", "reference": "CB", "axis_start": "CG", "axis_end": "SD", "blockers": {"CE"}, "bond_length": 1.79},
  ],
}

WATER_RESIDUE_NAMES = {"HOH", "WAT", "DOD", "OH2", "H2O"}
ODD_RESIDUE_CATEGORY_ORDER = (
  "Possible Misfits",
  "Unusual phi/psi in context",
  "Helical features",
  "Weak Sidechains",
  "Weak Backbone",
  "Weak regions",
  "Weak Waters/Ions",
  "Weak Ligands",
)


def _optional_import(module_name):
  try:
    return __import__(module_name)
  except Exception:
    return None


# Importing coot_gui eagerly at startup can trigger a second GTK/Rsvg stack in
# some bundled Coot builds (notably the CCP4 macOS bundle), which then breaks
# themed/icon rendering. Prefer the module Coot has already loaded, and only
# treat GUI Python helpers as available when that preloaded module exists.
def _coot_gui_module():
  existing_module = globals().get("coot_gui")
  if existing_module is not None:
    return existing_module
  existing_module = sys.modules.get("coot_gui")
  if existing_module is not None:
    globals()["coot_gui"] = existing_module
  return existing_module


coot_gui = _coot_gui_module()
coot_fitting = _optional_import("coot_fitting")
coot_ncs = _optional_import("coot_ncs")
GUI_PYTHON_AVAILABLE = coot_gui is not None
DIALOG_KEYBINDINGS_AVAILABLE = GUI_PYTHON_AVAILABLE

try:
  from gui_add_linked_cho import add_module_carbohydrate_gui
  CARBOHYDRATE_GUI_AVAILABLE = True
except Exception:
  CARBOHYDRATE_GUI_AVAILABLE = False


def _gui_unavailable_message(feature_name):
  return (
    feature_name
    + " is unavailable because the embedded Python environment cannot import the Gtk bindings needed by this helper."
  )


def _coot11_user_colour_warning():
  info_dialog(
    "Custom per-residue colouring is currently disabled in this Coot 1.1 build.\n\n"
    "The available Python user-defined colouring API stores persistent colour selections on the molecule "
    "and contaminates ordinary display modes such as Colour by Atom.\n\n"
    "This needs a separate overlay/additional-representation implementation rather than direct colouring of the base molecule."
  )
  return None


def _register_key_binding_if(condition, name, key, thunk):
  if condition:
    add_key_binding(name, key, thunk)


def _coot_state_file_path():
  xdg_state_home = os.environ.get("XDG_STATE_HOME")
  if xdg_state_home:
    state_home = xdg_state_home
  else:
    state_home = os.path.join(os.path.expanduser("~"), ".local", "state", "Coot")
  return os.path.join(state_home, "0-coot.state.py")


def _ensure_startup_state_file():
  state_file = _coot_state_file_path()
  if os.path.exists(state_file):
    return
  try:
    save_state()
  except Exception:
    print("Unable to create initial Coot state file:", state_file)
    traceback.print_exc()


def _apply_startup_settings():
  """Apply the top-of-file defaults once during startup."""
  set_symmetry_colour(*STARTUP_SYMMETRY_COLOUR)
  set_use_primary_mouse_button_for_view_rotation(STARTUP_USE_LEFT_MOUSE_FOR_ROTATION)
  _ensure_startup_state_file()
  set_scroll_by_wheel_mouse(STARTUP_SCROLL_BY_WHEEL_MOUSE)
  set_auto_fit_best_rotamer_clash_flag(STARTUP_AUTO_FIT_BEST_ROTAMER_CLASH)
  set_refine_max_residues(STARTUP_REFINE_MAX_RESIDUES)
  set_rotamer_search_mode(STARTUP_ROTAMER_SEARCH_MODE)
  set_map_sampling_rate(STARTUP_MAP_SAMPLING_RATE)
  allow_duplicate_sequence_numbers()
  set_add_terminal_residue_n_phi_psi_trials(STARTUP_ADD_TERMINAL_PHI_PSI_TRIALS)
  set_matrix(STARTUP_REFINEMENT_WEIGHT)
  set_smooth_scroll_flag(STARTUP_SMOOTH_SCROLL)
  set_active_map_drag_flag(STARTUP_ACTIVE_MAP_DRAG)
  set_show_environment_distances(STARTUP_SHOW_ENVIRONMENT_DISTANCES)
  set_show_environment_distances_bumps(STARTUP_SHOW_ENVIRONMENT_BUMPS)
  set_show_environment_distances_h_bonds(STARTUP_SHOW_ENVIRONMENT_H_BONDS)
  set_environment_distances_distance_limits(*STARTUP_ENVIRONMENT_DISTANCE_LIMITS)
  set_map_radius(STARTUP_MAP_RADIUS)
  set_map_radius_em(STARTUP_MAP_RADIUS_EM)
  set_default_bond_thickness(STARTUP_DEFAULT_BOND_THICKNESS)
  try:
    # Variable bond thickness is an on/off mode.  Do not feed a separate
    # "variable thickness" value to set_default_bond_thickness(): that API
    # changes the ordinary default bond thickness.
    set_use_variable_bond_thickness(STARTUP_USE_VARIABLE_BOND_THICKNESS)
  except NameError:
    print("Your coot is a bit old... consider upgrading...")
  set_default_temperature_factor_for_new_atoms(STARTUP_DEFAULT_NEW_ATOM_B)
  set_add_terminal_residue_do_post_refine(STARTUP_ADD_TERMINAL_POST_REFINE)
  set_terminal_residue_do_rigid_body_refine(STARTUP_TERMINAL_RIGID_BODY_REFINE)
  set_mutate_auto_fit_do_post_refine(STARTUP_MUTATE_AUTO_FIT_POST_REFINE)
  set_symmetry_size(STARTUP_SYMMETRY_SIZE)
  set_nomenclature_errors_on_read(STARTUP_NOMENCLATURE_ERRORS_MODE)
  set_reorienting_next_residue_mode(STARTUP_REORIENTING_NEXT_RESIDUE_MODE)
  if STARTUP_CURMUDGEON_MODE:
    _call_optional_coot_api("curmudgeon_mode")
  _call_optional_coot_api(
    "set_user_defined_rotation_centre_crosshairs_size_scale_factor",
    STARTUP_ROTATION_CENTRE_CROSSHAIRS_SCALE
  )
  _call_optional_coot_api(
    "set_rotation_centre_cross_hairs_colour",
    *STARTUP_ROTATION_CENTRE_CROSSHAIRS_COLOUR
  )


def _active_residue_or_status():
  """Return the active residue, or show a short status message if none exists."""
  residue = active_residue()
  if residue:
    return residue
  add_status_bar_text("No active residue")
  return None


def _scrollable_map_or_status():
  """Return the current scroll-wheel map, or show a short status message."""
  map_id = scroll_wheel_map()
  if map_id != -1 and map_id in map_molecule_list():
    return map_id
  map_id = imol_refinement_map()
  if map_id != -1 and map_id in map_molecule_list():
    try:
      set_scroll_wheel_map(map_id)
      set_scrollable_map(map_id)
    except Exception:
      pass
    return map_id
  for map_id in map_molecule_list():
    if map_is_displayed(map_id):
      try:
        set_scroll_wheel_map(map_id)
        set_scrollable_map(map_id)
      except Exception:
        pass
      return map_id
  map_ids = map_molecule_list()
  if map_ids:
    map_id = map_ids[0]
    try:
      set_scroll_wheel_map(map_id)
      set_scrollable_map(map_id)
    except Exception:
      pass
    return map_id
  add_status_bar_text("No active map")
  return None


def _active_model_molecule_for_map_tools(action_name):
  """Resolve the current model molecule for map/model wrapper helpers."""
  residue = active_residue()
  if residue and residue[0] in model_molecule_list():
    return residue[0]

  mol_id = go_to_atom_molecule_number()
  if mol_id in model_molecule_list():
    return mol_id

  add_status_bar_text(f"No active model for {action_name}")
  return None


def _active_chain_context_or_status():
  """Return the active residue's chain context as a stable snapshot."""
  residue = _active_residue_or_status()
  if not residue:
    return None
  return {
    "mol_id": residue[0],
    "chain_id": residue[1],
    "resno": residue[2],
    "ins_code": residue[3],
  }


def _active_atom_context_or_status():
  """Return the active atom context as a stable snapshot."""
  residue = _active_residue_or_status()
  if not residue:
    return None
  return {
    "mol_id": residue[0],
    "chain_id": residue[1],
    "resno": residue[2],
    "ins_code": residue[3],
    "atom_name": residue[4],
    "alt_conf": residue[5] if len(residue) > 5 else "",
  }


def _segment_specs_for_residue(mol_id, chain_id, resno, ins_code="", rows=None):
  rows = _segment_chain_snapshot(mol_id, chain_id) if rows is None else rows
  index = next((i for i, row in enumerate(rows)
                if (row["number"], row["ins"]) == (resno, ins_code or "")), None)
  if index is None or rows[index]["name"] not in _SEGMENT_POLYMER_NAMES:
    return []
  first = last = index
  while first and _segment_rows_connected(rows[first-1], rows[first]):
    first -= 1
  while last+1 < len(rows) and _segment_rows_connected(rows[last], rows[last+1]):
    last += 1
  return [(chain_id, row["number"], row["ins"]) for row in rows[first:last+1]]


def _segment_range_for_residue(mol_id, chain_id, resno, ins_code=""):
  """Return the contiguous polymer segment containing a residue number."""
  specs = _segment_specs_for_residue(mol_id, chain_id, resno, ins_code)
  if specs:
    return (specs[0][1], specs[-1][1])
  return None


def _active_segment_context_or_status():
  """Return active atom/segment context for segment-level helpers."""
  atom_context = _active_atom_context_or_status()
  if not atom_context:
    return None
  specs = _segment_specs_for_residue(
    atom_context["mol_id"],
    atom_context["chain_id"],
    atom_context["resno"],
    atom_context["ins_code"],
  )
  if not specs:
    add_status_bar_text("No active segment")
    return None
  atom_context["segment_start"] = specs[0][1]
  atom_context["segment_end"] = specs[-1][1]
  atom_context["segment_specs"] = specs
  return atom_context


def _refresh_extra_restraints_display(mol_id):
  """Redraw existing restraints without deleting them or changing visibility."""
  coot.graphics_draw()


def _parse_int_entry(value, error_message, minimum_value=None):
  """Parse an integer entry, showing a dialog on invalid input."""
  try:
    parsed_value = int(str(value).strip())
  except Exception:
    info_dialog(error_message)
    return None
  if minimum_value is not None and parsed_value < minimum_value:
    info_dialog(error_message)
    return None
  return parsed_value


def _renumber_residue_range_or_dialog(mol_id, ch_id, start_res, last_res, offset, failure_message=None):
  """Run residue renumbering and surface failures to the user."""
  if failure_message is None:
    failure_message = (
      "Failed to renumber residues.\n\n"
      "This usually means the requested renumbering would create overlapping sequence numbering."
    )
  try:
    status = renumber_residue_range(mol_id, ch_id, start_res, last_res, int(offset))
  except Exception:
    info_dialog(failure_message)
    return 0
  if not status:
    info_dialog(failure_message)
    return 0
  _mark_model_changed(mol_id)
  remembered = globals().get("NAVIGATION_LAST_RESIDUE")
  if remembered and remembered["mol_id"] == mol_id and remembered["chain_id"] == ch_id and start_res <= remembered["resno"] <= last_res:
    remembered["resno"] += int(offset)
  return status


def _mask_map_by_all_atoms(map_id, mol_id, invert_flag):
  """Mask a map using an explicit all-atoms selection when possible."""
  # Some Coot masking paths appear to treat the convenience molecule-wide call
  # a bit differently for hetero content. Prefer the explicit all-atoms
  # selection string first, then fall back to the broader helper.
  for atom_selection in ("/ /", "//"):
    try:
      masked_map_id = mask_map_by_atom_selection(map_id, mol_id, atom_selection, invert_flag)
    except Exception:
      masked_map_id = -1
    if masked_map_id != -1 and masked_map_id in map_molecule_list():
      return masked_map_id
  return mask_map_by_molecule(map_id, mol_id, invert_flag)


def _residue_is_polymer(mol_id, chain_id, resno, ins_code):
  if does_residue_exist_p(mol_id, chain_id, resno, ins_code) == 0:
    return False
  return residue_name(mol_id, chain_id, resno, ins_code) in POLYMER_RESIDUE_NAMES


def _residue_serial_number(mol_id, chain_id, resno, ins_code):
  """Return the exact serial number for a residue, including insertion code."""
  remembered = globals().get("NAVIGATION_LAST_RESIDUE")
  if remembered and (remembered["mol_id"], remembered["chain_id"], remembered["resno"], remembered["ins_code"]) == (mol_id, chain_id, resno, ins_code or ""):
    serial = remembered.get("serial_number")
    if serial is not None and 0 <= serial < chain_n_residues(chain_id, mol_id) and seqnum_from_serial_number(mol_id, chain_id, serial) == resno and insertion_code_from_serial_number(mol_id, chain_id, serial) == (ins_code or ""):
      return serial
  all_specs = getattr(coot, "all_residues_with_serial_numbers_py", None)
  if callable(all_specs):
    for entry in all_specs(mol_id) or []:
      parsed = _coot_residue_spec_from_spec(entry[1:])
      if parsed == (chain_id, resno, ins_code or ""):
        return int(entry[0])
    return -1
  for serial_number in range(chain_n_residues(chain_id, mol_id)):
    if seqnum_from_serial_number(mol_id, chain_id, serial_number) != resno:
      continue
    if insertion_code_from_serial_number(mol_id, chain_id, serial_number) == (ins_code or ""):
      return serial_number
  return -1


def _default_navigation_target_xyz(atom_xyz, residue_name_here):
  """Choose a stable fallback atom when orientation cannot be computed."""
  if not atom_xyz:
    return None
  if len(residue_name_here) == 1:
    preferred_atoms = ["C1'", "P", "C4'", "C3'"]
  else:
    preferred_atoms = ["CA", "CB", "N"]
  for preferred_atom in preferred_atoms:
    if preferred_atom in atom_xyz:
      return atom_xyz[preferred_atom]
  return next(iter(atom_xyz.values()))


def _first_polymer_residue_in_chain(mol_id, chain_id):
  for serial_number in range(chain_n_residues(chain_id, mol_id)):
    resno = seqnum_from_serial_number(mol_id, chain_id, serial_number)
    ins_code = insertion_code_from_serial_number(mol_id, chain_id, serial_number)
    if _residue_is_polymer(mol_id, chain_id, resno, ins_code):
      return (chain_id, resno, ins_code)
  return None


def _nearest_polymer_residue_to_rotation_centre(mol_id, chain_id=None, max_distance=10.0):
  """Find the nearest polymer residue centre near the current rotation centre."""
  centre = _rotation_centre_xyz()
  best_residue = None
  best_distance_sq = max_distance * max_distance
  for current_chain_id, resno, ins_code in _residue_specs_near_point(mol_id, centre, max_distance):
    if chain_id is not None and current_chain_id != chain_id:
      continue
    if not _residue_is_polymer(mol_id, current_chain_id, resno, ins_code):
      continue
    xyz = _trimmed_atom_xyz_map(mol_id, current_chain_id, resno, ins_code)
    anchor = _default_navigation_target_xyz(xyz, residue_name(mol_id, current_chain_id, resno, ins_code))
    if anchor is None:
      continue
    distance_sq = _distance_sq(centre, anchor)
    if distance_sq <= best_distance_sq:
      best_distance_sq = distance_sq
      best_residue = (current_chain_id, resno, ins_code)
  return best_residue


def _navigation_residue_near_rotation_centre(reference_residue, max_distance=6.0, distance_tolerance=0.25):
  """Prefer the polymer residue actually under the current rotation centre.

  After renumbering, Coot can leave active_residue()/stored navigation state
  pointing at the old residue number even though the view is still centred on
  the correct residue. Re-resolve the residue from geometry when the current
  view clearly favours a different residue in the same chain.
  """
  mol_id = reference_residue["mol_id"]
  chain_id = reference_residue["chain_id"]
  resno = reference_residue["resno"]
  ins_code = reference_residue["ins_code"]
  rotation_centre = _rotation_centre_xyz()
  reference_xyz = _trimmed_atom_xyz_map(mol_id, chain_id, resno, ins_code)
  reference_centre = _default_navigation_target_xyz(reference_xyz, residue_name(mol_id, chain_id, resno, ins_code))
  if reference_centre is not None and _distance_sq(rotation_centre, reference_centre) <= 0.5**2:
    return reference_residue
  nearest_residue = _nearest_polymer_residue_to_rotation_centre(mol_id, chain_id, max_distance=max_distance)
  if nearest_residue is None:
    return reference_residue
  if nearest_residue == (chain_id, resno, ins_code):
    return reference_residue
  if not isinstance(reference_centre, (list, tuple)) or len(reference_centre) != 3:
    reference_residue = dict(reference_residue)
    reference_residue.update({
      "chain_id": nearest_residue[0],
      "resno": nearest_residue[1],
      "ins_code": nearest_residue[2],
      "source": "nearby",
    })
    reference_residue.pop("serial_number", None)
    return reference_residue
  reference_distance_sq = _distance_sq(rotation_centre, reference_centre)
  nearest_xyz = _trimmed_atom_xyz_map(mol_id, *nearest_residue)
  nearest_centre = _default_navigation_target_xyz(nearest_xyz, residue_name(mol_id, *nearest_residue))
  if not isinstance(nearest_centre, (list, tuple)) or len(nearest_centre) != 3:
    return reference_residue
  nearest_distance_sq = _distance_sq(rotation_centre, nearest_centre)
  if nearest_distance_sq > 0.5**2:
    return reference_residue
  if nearest_distance_sq + (distance_tolerance * distance_tolerance) >= reference_distance_sq:
    return reference_residue
  reference_residue = dict(reference_residue)
  reference_residue.update({
    "chain_id": nearest_residue[0],
    "resno": nearest_residue[1],
    "ins_code": nearest_residue[2],
    "source": "nearby",
  })
  reference_residue.pop("serial_number", None)
  return reference_residue


def _navigation_reference_residue():
  """Resolve the residue that navigation should use as its starting point."""
  remembered = globals().get("NAVIGATION_LAST_RESIDUE")
  if remembered and remembered.get("view_centre") is not None and _distance_sq(_rotation_centre_xyz(), remembered["view_centre"]) < 1e-6:
    mol_id = remembered["mol_id"]
    if mol_id in model_molecule_list() and _residue_is_polymer(mol_id, remembered["chain_id"], remembered["resno"], remembered["ins_code"]):
      active = active_residue()
      spec = (remembered["chain_id"], remembered["resno"], remembered["ins_code"])
      anchor = _default_navigation_target_xyz(_trimmed_atom_xyz_map(mol_id, *spec), residue_name(mol_id, *spec))
      if anchor is not None and _distance_sq(anchor, remembered["view_centre"]) < 0.5**2 and (not active or active[0] == mol_id):
        result = dict(remembered)
        result["source"] = "stored"
        result["serial_number"] = _residue_serial_number(mol_id, remembered["chain_id"], remembered["resno"], remembered["ins_code"])
        return result
  residue = active_residue()
  if residue and _residue_is_polymer(residue[0], residue[1], residue[2], residue[3]):
    return _navigation_residue_near_rotation_centre({
      "mol_id": residue[0],
      "chain_id": residue[1],
      "resno": residue[2],
      "ins_code": residue[3],
      "source": "active",
    })

  global NAVIGATION_LAST_RESIDUE
  if NAVIGATION_LAST_RESIDUE:
    mol_id = NAVIGATION_LAST_RESIDUE["mol_id"]
    chain_id = NAVIGATION_LAST_RESIDUE["chain_id"]
    resno = NAVIGATION_LAST_RESIDUE["resno"]
    ins_code = NAVIGATION_LAST_RESIDUE["ins_code"]
    serial_number = NAVIGATION_LAST_RESIDUE.get("serial_number")
    if mol_id in model_molecule_list() and _residue_is_polymer(mol_id, chain_id, resno, ins_code):
      if serial_number is not None:
        if serial_number < 0 or serial_number >= chain_n_residues(chain_id, mol_id):
          serial_number = None
        elif seqnum_from_serial_number(mol_id, chain_id, serial_number) != resno:
          serial_number = None
        elif insertion_code_from_serial_number(mol_id, chain_id, serial_number) != (ins_code or ""):
          serial_number = None
      return _navigation_residue_near_rotation_centre({
        "mol_id": mol_id,
        "chain_id": chain_id,
        "resno": resno,
        "ins_code": ins_code,
        "source": "stored",
        "serial_number": serial_number,
      })
    NAVIGATION_LAST_RESIDUE = None

  mol_id = go_to_atom_molecule_number()
  model_molecules = model_molecule_list()
  if mol_id not in model_molecules:
    add_status_bar_text("No active polymer residue or model molecule")
    return None

  molecule_chain_ids = chain_ids(mol_id)
  chain_id = go_to_atom_chain_id()
  if chain_id not in molecule_chain_ids:
    chain_id = None

  nearby_residue = _nearest_polymer_residue_to_rotation_centre(mol_id, chain_id)
  if nearby_residue:
    return {
      "mol_id": mol_id,
      "chain_id": nearby_residue[0],
      "resno": nearby_residue[1],
      "ins_code": nearby_residue[2],
      "source": "nearby",
    }

  if chain_id:
    first_polymer = _first_polymer_residue_in_chain(mol_id, chain_id)
    if first_polymer:
      return {
        "mol_id": mol_id,
        "chain_id": first_polymer[0],
        "resno": first_polymer[1],
        "ins_code": first_polymer[2],
        "source": "chain_start",
      }

  for fallback_chain_id in molecule_chain_ids:
    first_polymer = _first_polymer_residue_in_chain(mol_id, fallback_chain_id)
    if first_polymer:
      return {
        "mol_id": mol_id,
        "chain_id": first_polymer[0],
        "resno": first_polymer[1],
        "ins_code": first_polymer[2],
        "source": "chain_start",
      }

  add_status_bar_text("No polymer residues found in the current model")
  return None


def _molecule_id_name_label(mol_id):
  molecule_label = molecule_name_stub_py(mol_id, 0) or molecule_name(mol_id) or ""
  molecule_label = os.path.basename(str(molecule_label).strip()) or f"molecule_{mol_id}"
  return f"Molecule #{mol_id}, {molecule_label}"


def _atom_b_factor_value(b_factor):
  if isinstance(b_factor, (list, tuple)) and b_factor:
    return _atom_b_factor_value(b_factor[0])
  try:
    return float(b_factor)
  except Exception:
    return None


def _residue_representative_atom_metrics(mol_id, chain_id, resno, ins_code):
  preferred_atom_names = ("CA", "P", "C1'", "C4'", "N", "C")
  atom_info = residue_info_py(mol_id, chain_id, resno, ins_code)
  if not atom_info:
    return None

  atom_records = []
  for atom in atom_info:
    try:
      atom_name = str(atom[0][0]).strip()
      alt_conf = str(atom[0][1]).strip()
      atom_details = atom[1]
      occupancy = float(atom_details[0])
      b_factor = _atom_b_factor_value(atom_details[1])
      element = str(atom_details[2]).strip().upper()
    except Exception:
      continue
    if not atom_name:
      continue
    atom_records.append(
      {
        "name": atom_name,
        "alt_conf": alt_conf,
        "occupancy": occupancy,
        "b_factor": b_factor,
        "element": element,
      }
    )

  if not atom_records:
    return None

  def choose_atom(predicate):
    for atom_record in atom_records:
      if predicate(atom_record):
        return atom_record
    return None

  representative_atom = None
  for preferred_atom_name in preferred_atom_names:
    representative_atom = choose_atom(
      lambda atom_record, preferred_atom_name=preferred_atom_name: (
        atom_record["name"] == preferred_atom_name and not atom_record["alt_conf"]
      )
    )
    if representative_atom:
      break
  if representative_atom is None:
    for preferred_atom_name in preferred_atom_names:
      representative_atom = choose_atom(
        lambda atom_record, preferred_atom_name=preferred_atom_name: atom_record["name"] == preferred_atom_name
      )
      if representative_atom:
        break
  if representative_atom is None:
    representative_atom = choose_atom(
      lambda atom_record: atom_record["element"] not in ("H", "D") and not atom_record["alt_conf"]
    )
  if representative_atom is None:
    representative_atom = choose_atom(lambda atom_record: atom_record["element"] not in ("H", "D"))
  if representative_atom is None:
    representative_atom = atom_records[0]
  return representative_atom


def _navigation_status_bar_label(mol_id, chain_id, resno, ins_code, residue_name_here=""):
  residue_label = _odd_residue_dialog_label(chain_id, resno, ins_code, residue_name_here)
  status_label = f"Current residue: {residue_label} in {_molecule_id_name_label(mol_id)}"
  representative_atom = _residue_representative_atom_metrics(mol_id, chain_id, resno, ins_code)
  if representative_atom is None:
    return status_label
  atom_name = representative_atom["name"]
  occupancy = representative_atom["occupancy"]
  b_factor = representative_atom["b_factor"]
  if b_factor is None:
    return f"{status_label} [{atom_name} occ {occupancy:.2f}]"
  return f"{status_label} [{atom_name} occ {occupancy:.2f}, B {b_factor:.1f}]"


def _go_to_navigation_residue(mol_id, chain_id, resno, ins_code, serial_number=None):
  """Navigate quietly by orienting/centring without invoking Go To highlighting."""
  global NAVIGATION_LAST_RESIDUE
  residue_name_here = residue_name(mol_id, chain_id, resno, ins_code)
  atom_xyz = _trimmed_atom_xyz_map(mol_id, chain_id, resno, ins_code)
  if serial_number is None:
    serial_number = _residue_serial_number(mol_id, chain_id, resno, ins_code)
  if not _orient_navigation_view(mol_id, chain_id, resno, ins_code, atom_xyz, residue_name_here):
    target_xyz = _default_navigation_target_xyz(atom_xyz, residue_name_here)
    if target_xyz is None:
      target_xyz = residue_centre_py(mol_id, chain_id, resno, ins_code)
    if isinstance(target_xyz, (list, tuple)) and len(target_xyz) == 3:
      set_rotation_centre(*target_xyz)
  NAVIGATION_LAST_RESIDUE = {
    "mol_id": mol_id,
    "chain_id": chain_id,
    "resno": resno,
    "ins_code": ins_code,
    "serial_number": serial_number,
    "view_centre": _rotation_centre_xyz(),
  }
  add_status_bar_text(
    _navigation_status_bar_label(mol_id, chain_id, resno, ins_code, residue_name_here)
  )
  return 1


def _rotation_centre_xyz():
  return [rotation_centre_position(axis) for axis in (0, 1, 2)]


def _distance_sq(point_1, point_2):
  dx = point_1[0] - point_2[0]
  dy = point_1[1] - point_2[1]
  dz = point_1[2] - point_2[2]
  return dx*dx + dy*dy + dz*dz


def _spatial_index(items, position, cell_size):
  bins = {}
  for item in items:
    xyz = position(item)
    key = tuple(math.floor(float(value) / cell_size) for value in xyz)
    bins.setdefault(key, []).append((xyz, item))
  return (float(cell_size), bins)


def _model_coordinate_snapshot(mol_id):
  """Read full-precision coordinates once; fall back on older Coot builds."""
  bulk = getattr(coot, "python_representation_kk", None)
  model_count = getattr(coot, "n_models", None)
  if callable(model_count) and model_count(mol_id) != 1:
    raise ValueError("This operation requires a single-model molecule.")
  if callable(bulk):
    models = bulk(mol_id)
    if not isinstance(models, (list, tuple)) or not models:
      raise RuntimeError("Coot could not read the model coordinates.")
    if len(models) != 1:
      raise ValueError("This operation requires a single-model molecule.")
    records = {}
    for chain, residues in models[0]:
      for number, ins, name, atoms in residues:
        spec = (str(chain), int(number), str(ins or ""))
        if spec in records:
          raise ValueError("Ambiguous residue identifiers in model snapshot.")
        records[spec] = {"name": name, "atoms": atoms}
    return records
  return {tuple(spec): {"name": residue_name(mol_id, *spec),
                        "atoms": residue_info_py(mol_id, *spec) or []}
          for spec in _all_residue_specs_for_colouring(mol_id)}


def _coordinate_snapshot_centres(snapshot):
  result = {}
  for spec, record in snapshot.items():
    points = [atom[2] for atom in record["atoms"]]
    if points:
      result[spec] = tuple(sum(point[axis] for point in points) / len(points) for axis in range(3))
  return result


def _spatial_query(index, point, radius):
  cell_size, bins = index
  lower = [math.floor((value-radius) / cell_size) for value in point]
  upper = [math.floor((value+radius) / cell_size) for value in point]
  radius_sq = radius * radius
  for x in range(lower[0], upper[0]+1):
    for y in range(lower[1], upper[1]+1):
      for z in range(lower[2], upper[2]+1):
        for xyz, item in bins.get((x, y, z), ()):
          if _distance_sq(xyz, point) <= radius_sq:
            yield item


def _residue_specs_near_point(mol_id, point, radius):
  native_query = getattr(coot, "residues_near_position_py", None)
  if callable(native_query):
    raw = native_query(mol_id, list(point), float(radius))
    if raw is False or raw is None:
      raise RuntimeError("Coot could not query nearby residues.")
    return list(dict.fromkeys(spec for spec in
      (_coot_residue_spec_from_spec(value) for value in raw) if spec is not None))
  specs = []
  for spec, record in _model_coordinate_snapshot(mol_id).items():
    if any(_distance_sq(atom[2], point) <= radius * radius
           for atom in record["atoms"]):
      specs.append(tuple(spec))
  return specs


def _vector_subtract(point_1, point_2):
  return [point_1[0] - point_2[0], point_1[1] - point_2[1], point_1[2] - point_2[2]]


def _vector_dot(vector_1, vector_2):
  return vector_1[0]*vector_2[0] + vector_1[1]*vector_2[1] + vector_1[2]*vector_2[2]


def _vector_cross(vector_1, vector_2):
  return [
    vector_1[1]*vector_2[2] - vector_1[2]*vector_2[1],
    vector_1[2]*vector_2[0] - vector_1[0]*vector_2[2],
    vector_1[0]*vector_2[1] - vector_1[1]*vector_2[0],
  ]


def _vector_length(vector):
  return math.sqrt(_vector_dot(vector, vector))


def _normalize_vector(vector):
  length = math.sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2])
  if length < 1.0e-7:
    return None
  return [component / length for component in vector]


def _scale_vector(vector, scale):
  return [component * scale for component in vector]


def _project_vector_perpendicular(vector, axis):
  axis_unit = _normalize_vector(axis)
  if axis_unit is None:
    return None
  projection_scale = _vector_dot(vector, axis_unit)
  perpendicular = [
    vector[0] - projection_scale * axis_unit[0],
    vector[1] - projection_scale * axis_unit[1],
    vector[2] - projection_scale * axis_unit[2],
  ]
  return _normalize_vector(perpendicular)


def _trimmed_atom_xyz_map(mol_id, chain_id, resno, ins_code):
  atom_info = residue_info_py(mol_id, chain_id, resno, ins_code)
  if not isinstance(atom_info, list):
    return {}
  atom_xyz = {}
  for atom in atom_info:
    try:
      atom_name = atom[0][0].strip()
      xyz = atom[2]
    except Exception:
      continue
    if atom_name and atom_name not in atom_xyz:
      atom_xyz[atom_name] = xyz
  return atom_xyz


def _rotation_matrix_to_quaternion(rotation_matrix):
  """Convert a 3x3 rotation matrix to Coot's view quaternion format."""
  m11, m12, m13 = rotation_matrix[0]
  m21, m22, m23 = rotation_matrix[1]
  m31, m32, m33 = rotation_matrix[2]
  trace = m11 + m22 + m33
  if trace > 0.0:
    scale = math.sqrt(trace + 1.0) * 2.0
    qw = 0.25 * scale
    qx = (m32 - m23) / scale
    qy = (m13 - m31) / scale
    qz = (m21 - m12) / scale
  elif m11 > m22 and m11 > m33:
    scale = math.sqrt(1.0 + m11 - m22 - m33) * 2.0
    qw = (m32 - m23) / scale
    qx = 0.25 * scale
    qy = (m12 + m21) / scale
    qz = (m13 + m31) / scale
  elif m22 > m33:
    scale = math.sqrt(1.0 + m22 - m11 - m33) * 2.0
    qw = (m13 - m31) / scale
    qx = (m12 + m21) / scale
    qy = 0.25 * scale
    qz = (m23 + m32) / scale
  else:
    scale = math.sqrt(1.0 + m33 - m11 - m22) * 2.0
    qw = (m21 - m12) / scale
    qx = (m13 + m31) / scale
    qy = (m23 + m32) / scale
    qz = 0.25 * scale
  return [qx, qy, qz, qw]


def _protein_navigation_orientation(atom_xyz):
  if "CA" not in atom_xyz:
    return None
  anchor = atom_xyz["CA"]
  if "CB" in atom_xyz and "N" in atom_xyz and "C" in atom_xyz:
    neighbour_vectors = []
    for atom_name in ("N", "C", "CB"):
      vector = _normalize_vector(_vector_subtract(atom_xyz[atom_name], anchor))
      if vector is not None:
        neighbour_vectors.append(vector)
    if len(neighbour_vectors) == 3:
      implicit_h = _normalize_vector([
        -(neighbour_vectors[0][0] + neighbour_vectors[1][0] + neighbour_vectors[2][0]),
        -(neighbour_vectors[0][1] + neighbour_vectors[1][1] + neighbour_vectors[2][1]),
        -(neighbour_vectors[0][2] + neighbour_vectors[1][2] + neighbour_vectors[2][2]),
      ])
      if implicit_h is not None:
        forward = _scale_vector(implicit_h, -1.0)
        up_guess = _project_vector_perpendicular(_vector_subtract(atom_xyz["CB"], anchor), forward)
        if up_guess is not None:
          return anchor, forward, up_guess
  if "N" in atom_xyz and "C" in atom_xyz and "O" in atom_xyz:
    neighbour_vectors = []
    for atom_name in ("N", "C", "O"):
      vector = _normalize_vector(_vector_subtract(atom_xyz[atom_name], anchor))
      if vector is not None:
        neighbour_vectors.append(vector)
    if len(neighbour_vectors) >= 2:
      pseudo_h = _normalize_vector([
        -sum(vector[0] for vector in neighbour_vectors),
        -sum(vector[1] for vector in neighbour_vectors),
        -sum(vector[2] for vector in neighbour_vectors),
      ])
      if pseudo_h is not None:
        forward = _scale_vector(pseudo_h, -1.0)
        up_guess = _project_vector_perpendicular(_vector_subtract(atom_xyz["O"], anchor), forward)
        if up_guess is not None:
          return anchor, forward, up_guess
  return None


def _glycine_navigation_orientation(atom_xyz):
  """Treat glycine like an alanine by constructing a consistent pseudo-CB direction."""
  if not all(atom_name in atom_xyz for atom_name in ("N", "CA", "C")):
    return None
  anchor = atom_xyz["CA"]
  n_vector = _normalize_vector(_vector_subtract(atom_xyz["N"], anchor))
  c_vector = _normalize_vector(_vector_subtract(atom_xyz["C"], anchor))
  if n_vector is None or c_vector is None:
    return None

  # Glycine has two alpha hydrogens. Construct the idealized tetrahedral pair
  # from the observed N-CA and C-CA directions, then treat one as a pseudo-CB
  # so glycine feels like the other amino acids during residue stepping.
  hydrogen_midpoint = [
    -(n_vector[0] + c_vector[0]) * 0.5,
    -(n_vector[1] + c_vector[1]) * 0.5,
    -(n_vector[2] + c_vector[2]) * 0.5,
  ]
  plane_normal = _normalize_vector(_vector_cross(n_vector, c_vector))
  if plane_normal is None:
    return None
  if "O" in atom_xyz:
    o_vector = _normalize_vector(_vector_subtract(atom_xyz["O"], anchor))
    if o_vector is not None and _vector_dot(plane_normal, o_vector) > 0.0:
      plane_normal = _scale_vector(plane_normal, -1.0)
  perpendicular_scale_sq = max(0.0, 1.0 - _vector_dot(hydrogen_midpoint, hydrogen_midpoint))
  perpendicular_scale = math.sqrt(perpendicular_scale_sq)
  pseudo_cb = _normalize_vector([
    hydrogen_midpoint[0] + perpendicular_scale * plane_normal[0],
    hydrogen_midpoint[1] + perpendicular_scale * plane_normal[1],
    hydrogen_midpoint[2] + perpendicular_scale * plane_normal[2],
  ])
  implicit_h = _normalize_vector([
    hydrogen_midpoint[0] - perpendicular_scale * plane_normal[0],
    hydrogen_midpoint[1] - perpendicular_scale * plane_normal[1],
    hydrogen_midpoint[2] - perpendicular_scale * plane_normal[2],
  ])
  if pseudo_cb is None or implicit_h is None:
    return None
  forward = _scale_vector(implicit_h, -1.0)
  up_guess = _project_vector_perpendicular(pseudo_cb, forward)
  if up_guess is None:
    return None
  return anchor, forward, up_guess


def _nucleotide_navigation_orientation(atom_xyz):
  if "C1'" not in atom_xyz:
    return None
  anchor = atom_xyz["C1'"]
  base_atom_name = None
  for atom_name in ("N9", "N1", "C4", "C2"):
    if atom_name in atom_xyz:
      base_atom_name = atom_name
      break
  if not base_atom_name:
    return None
  neighbour_vectors = []
  for atom_name in ("O4'", "C2'", base_atom_name):
    if atom_name not in atom_xyz:
      continue
    vector = _normalize_vector(_vector_subtract(atom_xyz[atom_name], anchor))
    if vector is not None:
      neighbour_vectors.append(vector)
  if len(neighbour_vectors) < 3:
    return None
  implicit_h = _normalize_vector([
    -sum(vector[0] for vector in neighbour_vectors),
    -sum(vector[1] for vector in neighbour_vectors),
    -sum(vector[2] for vector in neighbour_vectors),
  ])
  if implicit_h is None:
    return None
  forward = _scale_vector(implicit_h, -1.0)
  up_guess = _project_vector_perpendicular(_vector_subtract(atom_xyz[base_atom_name], anchor), forward)
  if up_guess is None:
    return None
  return anchor, forward, up_guess


def _nucleotide_base_plane_normal(atom_xyz):
  """Return a base-plane normal for canonical nucleotides when enough ring atoms are present."""
  base_plane_atom_sets = (
    ("N9", "C4", "C8"),  # purines
    ("N1", "C2", "C6"),  # pyrimidines
    ("C4", "C5", "C6"),  # generic aromatic fallback
  )
  for atom_names in base_plane_atom_sets:
    if not all(atom_name in atom_xyz for atom_name in atom_names):
      continue
    point_1 = atom_xyz[atom_names[0]]
    point_2 = atom_xyz[atom_names[1]]
    point_3 = atom_xyz[atom_names[2]]
    vector_1 = _normalize_vector(_vector_subtract(point_2, point_1))
    vector_2 = _normalize_vector(_vector_subtract(point_3, point_1))
    if vector_1 is None or vector_2 is None:
      continue
    plane_normal = _normalize_vector(_vector_cross(vector_1, vector_2))
    if plane_normal is not None:
      return plane_normal
  return None


def _orient_navigation_view(mol_id, chain_id, resno, ins_code, atom_xyz=None, residue_name_here=None):
  """Orient the residue for a sidechain/base-forward inspection view."""
  if atom_xyz is None:
    atom_xyz = _trimmed_atom_xyz_map(mol_id, chain_id, resno, ins_code)
  if not atom_xyz:
    return None
  if residue_name_here is None:
    residue_name_here = residue_name(mol_id, chain_id, resno, ins_code)
  if len(residue_name_here) == 1:
    orientation = _nucleotide_navigation_orientation(atom_xyz)
  elif residue_name_here == "GLY":
    orientation = _glycine_navigation_orientation(atom_xyz)
  else:
    orientation = _protein_navigation_orientation(atom_xyz)
  if not orientation:
    return None

  anchor, forward, up_guess = orientation
  right = _normalize_vector(_vector_cross(forward, up_guess))
  if right is None:
    return None
  up = _normalize_vector(_vector_cross(right, forward))
  if up is None:
    return None

  # For nucleotides, keep the vertical orientation but yaw so the aromatic base
  # sits more nearly in the screen plane rather than being obliquely foreshortened.
  if len(residue_name_here) == 1:
    base_plane_normal = _nucleotide_base_plane_normal(atom_xyz)
    if base_plane_normal is not None:
      desired_forward = _project_vector_perpendicular(base_plane_normal, up)
      if desired_forward is not None:
        if _vector_dot(desired_forward, forward) < 0.0:
          desired_forward = _scale_vector(desired_forward, -1.0)
        adjusted_right = _normalize_vector(_vector_cross(desired_forward, up))
        if adjusted_right is not None:
          adjusted_up = _normalize_vector(_vector_cross(adjusted_right, desired_forward))
          if adjusted_up is not None:
            forward = desired_forward
            right = adjusted_right
            up = adjusted_up

  rotation_matrix = [
    [right[0], right[1], right[2]],
    [up[0], up[1], up[2]],
    [-forward[0], -forward[1], -forward[2]],
  ]
  # Use the normal recenter path here so the local map recentres correctly.
  # The distracting green indicator turned out to come from the Go To sync,
  # which this custom navigation path no longer uses.
  set_rotation_centre(*anchor)
  set_view_quaternion(*_rotation_matrix_to_quaternion(rotation_matrix))
  return 1


def _capture_view_state():
  return {
    "rotation_centre": _rotation_centre_xyz(),
    "zoom": zoom_factor(),
    "quaternion": [get_view_quaternion_internal(i) for i in range(4)],
  }


def _restore_view_state(view_state):
  set_rotation_centre(*view_state["rotation_centre"])
  set_view_quaternion(*view_state["quaternion"])
  set_zoom(view_state["zoom"])


def _status_message(message):
  add_status_bar_text(message)


def _rotation_matrix_about_axis(axis, theta):
  normalized_axis = _normalize_vector(axis)
  if not normalized_axis:
    return None
  ux = normalized_axis[0]
  uy = normalized_axis[1]
  uz = normalized_axis[2]
  cos_theta = math.cos(theta)
  sin_theta = math.sin(theta)
  one_minus_cos = 1.0 - cos_theta
  return [
    cos_theta + ux * ux * one_minus_cos,
    ux * uy * one_minus_cos - uz * sin_theta,
    ux * uz * one_minus_cos + uy * sin_theta,
    uy * ux * one_minus_cos + uz * sin_theta,
    cos_theta + uy * uy * one_minus_cos,
    uy * uz * one_minus_cos - ux * sin_theta,
    uz * ux * one_minus_cos - uy * sin_theta,
    uz * uy * one_minus_cos + ux * sin_theta,
    cos_theta + uz * uz * one_minus_cos,
  ]


def _zeroify_rotation_matrix(matrix):
  if not matrix:
    return None
  zeroified = []
  for value in matrix:
    if abs(value) < 0.0000001:
      zeroified.append(0.0)
    elif abs(value - 1.0) < 0.0000001:
      zeroified.append(1.0)
    elif abs(value + 1.0) < 0.0000001:
      zeroified.append(-1.0)
    elif abs(value - 0.5) < 0.0000001:
      zeroified.append(0.5)
    elif abs(value + 0.5) < 0.0000001:
      zeroified.append(-0.5)
    else:
      zeroified.append(value)
  return zeroified


def _rotation_c2_xy(theta):
  return [
    math.cos(2.0 * theta), math.sin(2.0 * theta), 0.0,
    math.sin(2.0 * theta), -math.cos(2.0 * theta), 0.0,
    0.0, 0.0, -1.0,
  ]


def _icosahedral_phi():
  return 0.5 * (1.0 + math.sqrt(5.0))


def _icosahedral_axes_c2():
  phi = _icosahedral_phi()
  axes = []
  for sign_1 in [1.0, -1.0]:
    for sign_2 in [1.0, -1.0]:
      axes.append([0.0, sign_1, sign_2 * phi])
      axes.append([sign_1, sign_2 * phi, 0.0])
      axes.append([sign_2 * phi, 0.0, sign_1])
  return axes


def _icosahedral_axes_c3():
  axes = []
  for x in [1.0, -1.0]:
    for y in [1.0, -1.0]:
      for z in [1.0, -1.0]:
        axes.append([x, y, z])
  return axes


def _icosahedral_axes_c5():
  phi = _icosahedral_phi()
  axes = []
  for sign_1 in [1.0, -1.0]:
    for sign_2 in [1.0, -1.0]:
      axes.append([0.0, sign_1 / phi, sign_2 * phi])
      axes.append([sign_1 / phi, sign_2 * phi, 0.0])
      axes.append([sign_2 * phi, 0.0, sign_1 / phi])
  return axes


def _unique_rotation_matrices(matrices):
  unique = []
  seen = {}
  for matrix in matrices:
    zeroified = _zeroify_rotation_matrix(matrix)
    key = tuple([round(value, 6) for value in zeroified])
    if key not in seen:
      seen[key] = True
      unique.append(zeroified)
  return unique


def _generate_point_group_rotation_matrices(point_group_symbol):
  symbol = point_group_symbol.strip().upper()
  if not symbol:
    return None

  if symbol.startswith("C") and len(symbol) > 1 and symbol[1:].isdigit():
    order = int(symbol[1:])
    if order < 2:
      return None
    matrices = []
    for index in range(1, order):
      matrices.append(_rotation_matrix_about_axis([0.0, 0.0, 1.0], (2.0 * math.pi * index) / float(order)))
    return _unique_rotation_matrices(matrices)

  if symbol.startswith("D") and len(symbol) > 1 and symbol[1:].isdigit():
    order = int(symbol[1:])
    if order < 2:
      return None
    matrices = []
    for index in range(1, order):
      matrices.append(_rotation_matrix_about_axis([0.0, 0.0, 1.0], (2.0 * math.pi * index) / float(order)))
    for index in range(order):
      matrices.append(_rotation_c2_xy((math.pi * index) / float(order)))
    return _unique_rotation_matrices(matrices)

  if symbol == "T":
    axes = [
      [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],
      [1.0, 1.0, 1.0], [-1.0, -1.0, 1.0], [1.0, -1.0, -1.0], [-1.0, 1.0, -1.0],
    ]
    matrices = []
    for axis in axes[:3]:
      matrices.append(_rotation_matrix_about_axis(axis, math.pi))
    for axis in axes[3:]:
      matrices.append(_rotation_matrix_about_axis(axis, 2.0 * math.pi / 3.0))
      matrices.append(_rotation_matrix_about_axis(axis, 4.0 * math.pi / 3.0))
    return _unique_rotation_matrices(matrices)

  if symbol == "O":
    axes_c4 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    axes_c3 = [[1.0, 1.0, 1.0], [-1.0, -1.0, 1.0], [1.0, -1.0, -1.0], [-1.0, 1.0, -1.0]]
    axes_c2 = [
      [1.0, 1.0, 0.0], [1.0, -1.0, 0.0], [1.0, 0.0, 1.0],
      [1.0, 0.0, -1.0], [0.0, 1.0, 1.0], [0.0, 1.0, -1.0],
    ]
    matrices = []
    for axis in axes_c4:
      matrices.append(_rotation_matrix_about_axis(axis, math.pi / 2.0))
      matrices.append(_rotation_matrix_about_axis(axis, math.pi))
      matrices.append(_rotation_matrix_about_axis(axis, 3.0 * math.pi / 2.0))
    for axis in axes_c3:
      matrices.append(_rotation_matrix_about_axis(axis, 2.0 * math.pi / 3.0))
      matrices.append(_rotation_matrix_about_axis(axis, 4.0 * math.pi / 3.0))
    for axis in axes_c2:
      matrices.append(_rotation_matrix_about_axis(axis, math.pi))
    return _unique_rotation_matrices(matrices)

  if symbol == "I":
    matrices = []
    for axis in _icosahedral_axes_c2():
      matrices.append(_rotation_matrix_about_axis(axis, math.pi))
    for axis in _icosahedral_axes_c3():
      matrices.append(_rotation_matrix_about_axis(axis, 2.0 * math.pi / 3.0))
      matrices.append(_rotation_matrix_about_axis(axis, 4.0 * math.pi / 3.0))
    for axis in _icosahedral_axes_c5():
      matrices.append(_rotation_matrix_about_axis(axis, 2.0 * math.pi / 5.0))
      matrices.append(_rotation_matrix_about_axis(axis, 4.0 * math.pi / 5.0))
      matrices.append(_rotation_matrix_about_axis(axis, 6.0 * math.pi / 5.0))
      matrices.append(_rotation_matrix_about_axis(axis, 8.0 * math.pi / 5.0))
    return _unique_rotation_matrices(matrices)

  return None


def _close_smart_copy_template():
  global SMART_COPY_TEMPLATE_IMOL
  if SMART_COPY_TEMPLATE_IMOL in molecule_number_list():
    close_molecule(SMART_COPY_TEMPLATE_IMOL)
  SMART_COPY_TEMPLATE_IMOL = None


def _set_smart_copy_template(imol):
  global SMART_COPY_TEMPLATE_IMOL
  SMART_COPY_TEMPLATE_IMOL = imol
  set_mol_displayed(imol, 0)
  set_mol_active(imol, 0)


def _smart_copy_atom_selection(chain_id, resno, ins_code):
  if ins_code:
    return None
  return "//{chain_id}/{resno}".format(chain_id=chain_id, resno=resno)


def _find_model_molecule_for_click_spec(chain_id, resno, ins_code):
  """Find the clicked model by residue identity when the click payload lacks a usable imol."""
  if not chain_id or resno is False:
    return -1
  matches = [imol for imol in model_molecule_list()
             if residue_exists_qm(imol, chain_id, resno, ins_code)]
  return matches[0] if len(matches) == 1 else -1


def _click_spec_field(click_spec, long_index, short_index, default=None):
  """Read a field from either of the common Coot 1.x click-spec layouts."""
  if not isinstance(click_spec, (list, tuple)):
    return default
  if len(click_spec) >= 7:
    return click_spec[long_index]
  if len(click_spec) == 6:
    return click_spec[short_index]
  return default


def _click_spec_imol(click_spec):
  """Decode molecule identity without confusing PDB model number with imol."""
  valid_model_mols = model_molecule_list()
  if not isinstance(click_spec, (list, tuple)):
    return -1
  imol = _coot_int_arg(_click_spec_field(click_spec, 1, 0, -1))
  if isinstance(imol, int) and imol in valid_model_mols:
    return imol
  if len(click_spec) >= 6:
    # An explicit but invalid imol must never redirect an edit to another model.
    return -1
  chain_id = _click_spec_chain_id(click_spec)
  resno = _click_spec_res_no(click_spec)
  ins_code = _click_spec_ins_code(click_spec)
  imol = _find_model_molecule_for_click_spec(chain_id, resno, ins_code)
  if imol in valid_model_mols:
    return imol
  return -1


def _click_spec_chain_id(click_spec):
  chain_id = _click_spec_field(click_spec, 2, 1, False)
  return False if chain_id is False else _coot_text_arg(chain_id)


def _click_spec_res_no(click_spec):
  resno = _coot_int_arg(_click_spec_field(click_spec, 3, 2, False))
  return resno if resno is not None else False


def _click_spec_ins_code(click_spec):
  return _coot_text_arg(_click_spec_field(click_spec, 4, 3, ""))


def _click_spec_atom_name(click_spec):
  return _coot_text_arg(_click_spec_field(click_spec, 5, 4, ""))


def _click_spec_alt_conf(click_spec):
  return _coot_text_arg(_click_spec_field(click_spec, 6, 5, ""))


def _resolve_script_click(spec):
  """Repair GTK4's unset imol without guessing from the PDB model number or cofr."""
  if not isinstance(spec, (list, tuple)) or len(spec) not in (6, 7):
    raise ValueError("Coot returned an unrecognised atom pick; please click again.")
  index = 1 if len(spec) == 7 else 0
  imol = _coot_int_arg(spec[index])
  models = model_molecule_list()
  chain, number, ins = _click_spec_chain_id(spec), _click_spec_res_no(spec), _click_spec_ins_code(spec)
  name, alt = _click_spec_atom_name(spec).strip(), _click_spec_alt_conf(spec).strip()
  if chain is False or number is False or not name:
    raise ValueError("Coot returned an incomplete atom pick; please click again.")
  if imol != -1 and imol not in models:
    raise ValueError("The clicked model is no longer available; please click again.")
  displayed = getattr(coot, "mol_is_displayed", None)
  pickable = getattr(coot, "mol_is_active", None)
  matches = []
  for model in models if imol == -1 else [imol]:
    if callable(displayed) and not displayed(model):
      continue
    if callable(pickable) and not pickable(model):
      continue
    if not residue_exists_qm(model, chain, number, ins):
      continue
    atoms = coot.residue_info_py(model, chain, number, ins) or []
    if any(a[0][0].strip() == name and str(a[0][1] or "").strip() == alt for a in atoms):
      matches.append(model)
  if len(matches) != 1:
    if matches:
      raise ValueError("Ambiguous atom pick: hide the other matching models, then click again.")
    raise ValueError("The clicked atom is not in a displayed, pickable model; please click again.")
  repaired = list(spec)
  repaired[index] = matches[0]
  return repaired


_SCRIPT_CLICK_STATE = {}


def _cancel_script_clicks():
  state = dict(_SCRIPT_CLICK_STATE)
  _SCRIPT_CLICK_STATE.clear()
  if state.get("picking"):
    clear = getattr(coot, "clear_pending_picks", None)
    if callable(clear):
      try:
        clear()
      except Exception as error:
        print("Could not cancel atom picking:", error)
  window = state.get("window")
  if window is not None:
    window.destroy()
  cleanup = state.get("cleanup")
  if callable(cleanup):
    cleanup()


def clear_pending_custom_picks():
  """Cancel script pickers too, so Tab cannot leave endpoint markers behind."""
  _cancel_script_clicks()
  _clear_annotation_range_pick()
  coot.clear_pending_picks()


def _script_user_defined_click(count, callback, on_pick=None, cleanup=None, title=None, prompt=None):
  """Resolve picks immediately, keeping picking instructions in the status bar."""
  _cancel_script_clicks()
  clear_annotation = globals().get("_clear_annotation_range_pick")
  if callable(clear_annotation):
    clear_annotation()
  token = object()
  picks = []
  _SCRIPT_CLICK_STATE.update(token=token, picking=True, cleanup=cleanup)
  def current():
    return _SCRIPT_CLICK_STATE.get("token") is token
  def cancel():
    if current():
      _cancel_script_clicks()
  def arm():
    coot.user_defined_click_py(1, picked)
  def picked(raw):
    if not current():
      return
    try:
      spec = _resolve_script_click(raw)
    except ValueError as error:
      add_status_bar_text(str(error))
      print("Atom picking:", error)
      try:
        arm()
      except Exception:
        cancel()
        raise
      return
    except Exception:
      cancel()
      raise
    try:
      picks.append(spec)
      if callable(on_pick):
        on_pick(spec, len(picks))
      if len(picks) < count:
        arm()
        return
      # A model may have been closed between the two clicks.
      if any(_click_spec_imol(pick) < 0 for pick in picks):
        raise ValueError("A picked model was closed; select the range again.")
      _SCRIPT_CLICK_STATE["picking"] = False
      window = _SCRIPT_CLICK_STATE.pop("window", None)
      if window is not None:
        window.destroy()
      return callback(*picks)
    except Exception:
      cancel()
      raise
  try:
    if prompt or title:
      add_status_bar_text(prompt or f"{title}. Tab cancels.")
    arm()
  except Exception:
    cancel()
    raise
  return token


def ncs_master_chain_id(imol):
  if hasattr(coot, "ncs_master_chains_py"):
    master_chains = coot.ncs_master_chains_py(imol)
    if isinstance(master_chains, list) and master_chains:
      first_master = master_chains[0]
      if isinstance(first_master, str) and first_master:
        return first_master
      if isinstance(first_master, (list, tuple)) and first_master:
        return first_master[0]
  ncs_groups = ncs_chain_ids(imol)
  if isinstance(ncs_groups, list) and ncs_groups:
    first_group = ncs_groups[0]
    if isinstance(first_group, list) and first_group:
      return first_group[0]
  return False


def refine_residues(imol, residue_specs):
  if not hasattr(coot, "refine_residues_py"):
    raise NameError("refine_residues_py is not available in this Coot build")
  previous_immediate_replacement = 0
  if hasattr(coot, "refinement_immediate_replacement_state"):
    previous_immediate_replacement = refinement_immediate_replacement_state()
  try:
    set_refinement_immediate_replacement(1)
    refine_residues_py(imol, residue_specs)
    return accept_moving_atoms_py()
  finally:
    set_refinement_immediate_replacement(previous_immediate_replacement)


if coot_fitting and hasattr(coot_fitting, "pepflip_active_residue"):
  pepflip_active_residue = coot_fitting.pepflip_active_residue


_preloaded_coot_gui = _coot_gui_module()

if _preloaded_coot_gui is not None:
  coot_toolbar_button = _preloaded_coot_gui.coot_toolbar_button
  attach_module_menu_button = _preloaded_coot_gui.attach_module_menu_button
  generic_button_dialog = _preloaded_coot_gui.generic_button_dialog
  generic_multiple_entries_with_check_button = (
    _preloaded_coot_gui.generic_multiple_entries_with_check_button
  )
  add_simple_action_to_menu = _preloaded_coot_gui.add_simple_action_to_menu
  add_module_cryo_em_gui = _preloaded_coot_gui.add_module_cryo_em_gui
  add_module_refine = _preloaded_coot_gui.add_module_refine
  Gtk = _preloaded_coot_gui.Gtk

  _menu_action_counter = 0

  def add_simple_coot_menu_menuitem(menu, menu_item_label, activate_function):
    if menu is None:
      return None
    global _menu_action_counter
    action_name = "coot_trimmings_action_{0}".format(_menu_action_counter)
    _menu_action_counter += 1

    def on_activate(_action, _parameter):
      return activate_function(None)

    add_simple_action_to_menu(menu, menu_item_label, action_name, on_activate)

  def generic_single_entry(function_label, entry_1_default_text, go_button_label, handle_go_function):
    """Small GTK4-safe fallback for the common one-text-entry prompt."""
    window = Gtk.Window()
    window.set_title("Coot")

    vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=8)
    hbox_buttons = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
    label = Gtk.Label(label=function_label)
    entry = Gtk.Entry()
    cancel_button = Gtk.Button(label="Cancel")
    go_button = Gtk.Button(label=go_button_label)

    label.set_margin_start(12)
    label.set_margin_end(12)
    label.set_margin_top(12)
    label.set_margin_bottom(4)
    entry.set_margin_start(12)
    entry.set_margin_end(12)
    entry.set_margin_bottom(8)
    hbox_buttons.set_margin_start(12)
    hbox_buttons.set_margin_end(12)
    hbox_buttons.set_margin_bottom(12)

    if isinstance(entry_1_default_text, str):
      entry.set_text(entry_1_default_text)

    def close_window(*_args):
      window.destroy()
      return False

    def submit(*_args):
      handle_go_function(entry.get_text())
      window.destroy()
      return False

    cancel_button.connect("clicked", close_window)
    go_button.connect("clicked", submit)
    entry.connect("activate", submit)

    vbox.append(label)
    vbox.append(entry)
    hbox_buttons.append(cancel_button)
    hbox_buttons.append(go_button)
    vbox.append(hbox_buttons)
    window.set_child(vbox)
    window.present()
    entry.grab_focus()

  def generic_confirm_dialog(
    function_label,
    message_text,
    cancel_button_label,
    handle_cancel_function,
    confirm_button_label,
    handle_confirm_function,
    handle_close_function=None,
  ):
    """Small GTK4-safe two-option dialog for potentially expensive actions."""
    window = Gtk.Window()
    window.set_title("Coot")
    try:
      window.set_modal(True)
    except Exception:
      pass

    vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=8)
    hbox_buttons = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
    title_label = Gtk.Label(label=function_label)
    label = Gtk.Label(label=message_text)
    cancel_button = Gtk.Button(label=cancel_button_label)
    confirm_button = Gtk.Button(label=confirm_button_label)

    title_label.set_margin_start(12)
    title_label.set_margin_end(12)
    title_label.set_margin_top(12)
    title_label.set_margin_bottom(4)
    title_label.set_wrap(True)
    title_label.set_xalign(0.0)
    label.set_margin_start(12)
    label.set_margin_end(12)
    label.set_margin_top(0)
    label.set_margin_bottom(8)
    label.set_wrap(True)
    label.set_xalign(0.0)
    hbox_buttons.set_margin_start(12)
    hbox_buttons.set_margin_end(12)
    hbox_buttons.set_margin_bottom(12)
    hbox_buttons.set_halign(Gtk.Align.CENTER)

    def close_window(*_args):
      window.destroy()
      return False

    def cancel(*_args):
      handle_cancel_function()
      window.destroy()
      return False

    def submit(*_args):
      handle_confirm_function()
      window.destroy()
      return False

    cancel_button.connect("clicked", cancel)
    confirm_button.connect("clicked", submit)
    if handle_close_function is not None:
      def on_close(*_args):
        handle_close_function()
        return False
      window.connect("close-request", on_close)

    vbox.append(title_label)
    vbox.append(label)
    hbox_buttons.append(cancel_button)
    hbox_buttons.append(confirm_button)
    vbox.append(hbox_buttons)
    window.set_child(vbox)
    window.present()

  def generic_double_entry(
    function_label,
    entry_1_label,
    entry_1_default_text,
    entry_2_label,
    entry_2_default_text,
    go_button_label,
    handle_go_function,
  ):
    """GTK4-safe two-entry prompt used for small add/edit helpers."""
    window = Gtk.Window()
    window.set_title("Coot")

    vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=8)
    hbox_buttons = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
    title_label = Gtk.Label(label=function_label)
    label_1 = Gtk.Label(label=entry_1_label)
    entry_1 = Gtk.Entry()
    label_2 = Gtk.Label(label=entry_2_label)
    entry_2 = Gtk.Entry()
    cancel_button = Gtk.Button(label="Cancel")
    go_button = Gtk.Button(label=go_button_label)

    title_label.set_margin_start(12)
    title_label.set_margin_end(12)
    title_label.set_margin_top(12)
    title_label.set_margin_bottom(4)
    label_1.set_margin_start(12)
    label_1.set_margin_end(12)
    entry_1.set_margin_start(12)
    entry_1.set_margin_end(12)
    label_2.set_margin_start(12)
    label_2.set_margin_end(12)
    label_2.set_margin_top(6)
    entry_2.set_margin_start(12)
    entry_2.set_margin_end(12)
    entry_2.set_margin_bottom(8)
    hbox_buttons.set_margin_start(12)
    hbox_buttons.set_margin_end(12)
    hbox_buttons.set_margin_bottom(12)

    if isinstance(entry_1_default_text, str):
      entry_1.set_text(entry_1_default_text)
    if isinstance(entry_2_default_text, str):
      entry_2.set_text(entry_2_default_text)

    def close_window(*_args):
      window.destroy()
      return False

    def submit(*_args):
      status = handle_go_function(entry_1.get_text(), entry_2.get_text())
      if status not in (0, False):
        window.destroy()
      return False

    cancel_button.connect("clicked", close_window)
    go_button.connect("clicked", submit)
    entry_1.connect("activate", lambda *_args: entry_2.grab_focus())
    entry_2.connect("activate", submit)

    vbox.append(title_label)
    vbox.append(label_1)
    vbox.append(entry_1)
    vbox.append(label_2)
    vbox.append(entry_2)
    hbox_buttons.append(cancel_button)
    hbox_buttons.append(go_button)
    vbox.append(hbox_buttons)
    window.set_child(vbox)
    window.present()
    entry_1.grab_focus()

  def _populate_entry_from_file_chooser(parent_window, target_entry, chooser_title,
                                      initial_directory=None, all_files=False):
    """Open a GTK file chooser and copy the selected path into an entry."""
    def _gtk_gio_module():
      gio_module = globals().get("Gio")
      if gio_module is None:
        gio_module = _coot_gui_repository_module("Gio")
        if gio_module is not None:
          globals()["Gio"] = gio_module
      if gio_module is None:
        return None
      globals()["Gio"] = gio_module
      return gio_module

    def _remember_browsed_directory_from_path(path):
      global COMMON_MONOMER_LAST_BROWSED_DIRECTORY
      if not isinstance(path, str):
        return None
      normalised_path = os.path.abspath(os.path.expanduser(path))
      directory = normalised_path if os.path.isdir(normalised_path) else os.path.dirname(normalised_path)
      if directory and os.path.isdir(directory):
        COMMON_MONOMER_LAST_BROWSED_DIRECTORY = directory
      return None

    def _apply_file_chooser_initial_folder(chooser, gio_module=None):
      directory = initial_directory or COMMON_MONOMER_LAST_BROWSED_DIRECTORY
      if not directory or not os.path.isdir(directory):
        return None
      if gio_module is None:
        gio_module = _gtk_gio_module()
      if gio_module is None or not hasattr(gio_module, "File"):
        return None
      initial_folder = gio_module.File.new_for_path(directory)
      if hasattr(chooser, "set_initial_folder"):
        try:
          chooser.set_initial_folder(initial_folder)
          return None
        except Exception:
          pass
      if hasattr(chooser, "set_current_folder"):
        try:
          with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            chooser.set_current_folder(initial_folder)
        except Exception:
          return None
      return None

    cif_filter = Gtk.FileFilter()
    cif_filter.set_name("CIF dictionaries")
    for pattern in ("*.cif", "*.mmcif", "*.dic"):
      cif_filter.add_pattern(pattern)

    all_filter = Gtk.FileFilter()
    all_filter.set_name("All files")
    all_filter.add_pattern("*")

    if hasattr(Gtk, "FileDialog"):
      gio_module = _gtk_gio_module()

      chooser = Gtk.FileDialog()
      chooser.set_title(chooser_title)
      _apply_file_chooser_initial_folder(chooser, gio_module)

      if gio_module is not None and hasattr(gio_module, "ListStore"):
        filters = gio_module.ListStore.new(Gtk.FileFilter)
        filters.append(cif_filter)
        filters.append(all_filter)
        chooser.set_filters(filters)
        chooser.set_default_filter(all_filter if all_files else cif_filter)

      def on_open_finished(dialog, result):
        try:
          chosen_file = dialog.open_finish(result)
        except Exception:
          return None
        if chosen_file is not None and chosen_file.get_path():
          chosen_path = chosen_file.get_path()
          target_entry.set_text(chosen_path)
          _remember_browsed_directory_from_path(chosen_path)
        return None

      chooser.open(parent_window, None, on_open_finished)
      return None

    if hasattr(Gtk, "FileChooserNative"):
      # Coot's embedded GTK can lack Gtk.FileDialog even though the old
      # chooser still works. Keep the fallback, but suppress the Python-side
      # deprecation warnings so the console does not fill up on each browse.
      with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        chooser = Gtk.FileChooserNative.new(
          chooser_title,
          parent_window,
          Gtk.FileChooserAction.OPEN,
          "Open",
          "Cancel",
        )
        chooser.add_filter(cif_filter)
        chooser.add_filter(all_filter)
        chooser.set_filter(all_filter if all_files else cif_filter)
      _apply_file_chooser_initial_folder(chooser)

      def on_response(dialog, response_id):
        if response_id == Gtk.ResponseType.ACCEPT:
          chosen_file = dialog.get_file()
          if chosen_file is not None and chosen_file.get_path():
            chosen_path = chosen_file.get_path()
            target_entry.set_text(chosen_path)
            _remember_browsed_directory_from_path(chosen_path)
        dialog.destroy()

      chooser.connect("response", on_response)
      chooser.show()
      return None

    if hasattr(Gtk, "FileChooserDialog"):
      with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        chooser = Gtk.FileChooserDialog(
          title=chooser_title,
          transient_for=parent_window,
          action=Gtk.FileChooserAction.OPEN,
        )
        chooser.add_button("Cancel", Gtk.ResponseType.CANCEL)
        chooser.add_button("Open", Gtk.ResponseType.ACCEPT)
        chooser.add_filter(cif_filter)
        chooser.add_filter(all_filter)
        chooser.set_filter(all_filter if all_files else cif_filter)
      _apply_file_chooser_initial_folder(chooser)

      def on_response(dialog, response_id):
        if response_id == Gtk.ResponseType.ACCEPT:
          chosen_file = dialog.get_file()
          if chosen_file is not None and chosen_file.get_path():
            chosen_path = chosen_file.get_path()
            target_entry.set_text(chosen_path)
            _remember_browsed_directory_from_path(chosen_path)
        dialog.destroy()

      chooser.connect("response", on_response)
      chooser.present()
      return None

    info_dialog("Browse is unavailable in this Gtk build; please paste the file path manually.")
    return None

  def _populate_entry_from_save_file_chooser(parent_window, target_entry, chooser_title):
    """Open a GTK save chooser and copy the selected output path into an entry."""
    def _gtk_gio_module():
      gio_module = globals().get("Gio")
      if gio_module is None:
        gio_module = _coot_gui_repository_module("Gio")
        if gio_module is not None:
          globals()["Gio"] = gio_module
      if gio_module is None:
        return None
      globals()["Gio"] = gio_module
      return gio_module

    current_path = os.path.abspath(os.path.expanduser(target_entry.get_text().strip() or ""))
    current_directory = os.path.dirname(current_path) if current_path else ""
    current_filename = os.path.basename(current_path) if current_path else ""

    cif_filter = Gtk.FileFilter()
    cif_filter.set_name("mmCIF files")
    for pattern in ("*.cif", "*.mmcif"):
      cif_filter.add_pattern(pattern)

    all_filter = Gtk.FileFilter()
    all_filter.set_name("All files")
    all_filter.add_pattern("*")

    if hasattr(Gtk, "FileDialog"):
      gio_module = _gtk_gio_module()
      chooser = Gtk.FileDialog()
      chooser.set_title(chooser_title)

      if gio_module is not None and hasattr(gio_module, "File") and current_directory and os.path.isdir(current_directory):
        try:
          chooser.set_initial_folder(gio_module.File.new_for_path(current_directory))
        except Exception:
          pass
      if current_filename and hasattr(chooser, "set_initial_name"):
        try:
          chooser.set_initial_name(current_filename)
        except Exception:
          pass

      if gio_module is not None and hasattr(gio_module, "ListStore"):
        filters = gio_module.ListStore.new(Gtk.FileFilter)
        filters.append(cif_filter)
        filters.append(all_filter)
        chooser.set_filters(filters)
        chooser.set_default_filter(cif_filter)

      def on_save_finished(dialog, result):
        try:
          chosen_file = dialog.save_finish(result)
        except Exception:
          return None
        if chosen_file is not None and chosen_file.get_path():
          target_entry.set_text(chosen_file.get_path())
        return None

      chooser.save(parent_window, None, on_save_finished)
      return None

    if hasattr(Gtk, "FileChooserNative"):
      with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        chooser = Gtk.FileChooserNative.new(
          chooser_title,
          parent_window,
          Gtk.FileChooserAction.SAVE,
          "Save",
          "Cancel",
        )
        chooser.add_filter(cif_filter)
        chooser.add_filter(all_filter)
        if hasattr(chooser, "set_do_overwrite_confirmation"):
          chooser.set_do_overwrite_confirmation(True)
        if current_directory and hasattr(chooser, "set_current_folder"):
          gio_module = _gtk_gio_module()
          if gio_module is not None and hasattr(gio_module, "File") and os.path.isdir(current_directory):
            try:
              chooser.set_current_folder(gio_module.File.new_for_path(current_directory))
            except Exception:
              pass
        if current_filename and hasattr(chooser, "set_current_name"):
          try:
            chooser.set_current_name(current_filename)
          except Exception:
            pass

      def on_response(dialog, response_id):
        if response_id == Gtk.ResponseType.ACCEPT:
          chosen_file = dialog.get_file()
          if chosen_file is not None and chosen_file.get_path():
            target_entry.set_text(chosen_file.get_path())
        dialog.destroy()

      chooser.connect("response", on_response)
      chooser.show()
      return None

    if hasattr(Gtk, "FileChooserDialog"):
      with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        chooser = Gtk.FileChooserDialog(
          title=chooser_title,
          transient_for=parent_window,
          action=Gtk.FileChooserAction.SAVE,
        )
        chooser.add_button("Cancel", Gtk.ResponseType.CANCEL)
        chooser.add_button("Save", Gtk.ResponseType.ACCEPT)
        chooser.add_filter(cif_filter)
        chooser.add_filter(all_filter)
        if hasattr(chooser, "set_do_overwrite_confirmation"):
          chooser.set_do_overwrite_confirmation(True)
        if current_directory and hasattr(chooser, "set_current_folder"):
          gio_module = _gtk_gio_module()
          if gio_module is not None and hasattr(gio_module, "File") and os.path.isdir(current_directory):
            try:
              chooser.set_current_folder(gio_module.File.new_for_path(current_directory))
            except Exception:
              pass
        if current_filename and hasattr(chooser, "set_current_name"):
          try:
            chooser.set_current_name(current_filename)
          except Exception:
            pass

      def on_response(dialog, response_id):
        if response_id == Gtk.ResponseType.ACCEPT:
          chosen_file = dialog.get_file()
          if chosen_file is not None and chosen_file.get_path():
            target_entry.set_text(chosen_file.get_path())
        dialog.destroy()

      chooser.connect("response", on_response)
      chooser.present()
      return None

    info_dialog("Browse is unavailable in this Gtk build; please paste the output path manually.")
    return None

  def generic_double_entry_with_file_browse(
    function_label,
    entry_1_label,
    entry_1_default_text,
    entry_2_label,
    entry_2_default_text,
    entry_2_browse_title,
    go_button_label,
    handle_go_function,
  ):
    """GTK4-safe two-entry prompt with a Browse button for the second field."""
    window = Gtk.Window()
    window.set_title("Coot")

    vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=8)
    hbox_buttons = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
    entry_2_row = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
    title_label = Gtk.Label(label=function_label)
    label_1 = Gtk.Label(label=entry_1_label)
    entry_1 = Gtk.Entry()
    label_2 = Gtk.Label(label=entry_2_label)
    entry_2 = Gtk.Entry()
    browse_button = Gtk.Button(label="Browse...")
    cancel_button = Gtk.Button(label="Cancel")
    go_button = Gtk.Button(label=go_button_label)

    title_label.set_margin_start(12)
    title_label.set_margin_end(12)
    title_label.set_margin_top(12)
    title_label.set_margin_bottom(4)
    label_1.set_margin_start(12)
    label_1.set_margin_end(12)
    entry_1.set_margin_start(12)
    entry_1.set_margin_end(12)
    label_2.set_margin_start(12)
    label_2.set_margin_end(12)
    label_2.set_margin_top(6)
    entry_2_row.set_margin_start(12)
    entry_2_row.set_margin_end(12)
    entry_2_row.set_margin_bottom(8)
    entry_2.set_hexpand(True)
    hbox_buttons.set_margin_start(12)
    hbox_buttons.set_margin_end(12)
    hbox_buttons.set_margin_bottom(12)

    if isinstance(entry_1_default_text, str):
      entry_1.set_text(entry_1_default_text)
    if isinstance(entry_2_default_text, str):
      entry_2.set_text(entry_2_default_text)

    def close_window(*_args):
      window.destroy()
      return False

    def submit(*_args):
      status = handle_go_function(entry_1.get_text(), entry_2.get_text())
      if status not in (0, False):
        window.destroy()
      return False

    cancel_button.connect("clicked", close_window)
    browse_button.connect(
      "clicked",
      lambda *_args: _populate_entry_from_file_chooser(window, entry_2, entry_2_browse_title),
    )
    go_button.connect("clicked", submit)
    entry_1.connect("activate", lambda *_args: entry_2.grab_focus())
    entry_2.connect("activate", submit)

    entry_2_row.append(entry_2)
    entry_2_row.append(browse_button)

    vbox.append(title_label)
    vbox.append(label_1)
    vbox.append(entry_1)
    vbox.append(label_2)
    vbox.append(entry_2_row)
    hbox_buttons.append(cancel_button)
    hbox_buttons.append(go_button)
    vbox.append(hbox_buttons)
    window.set_child(vbox)
    window.present()
    entry_1.grab_focus()

  def generic_single_entry_with_file_browse(
    function_label,
    entry_1_default_text,
    entry_1_browse_title,
    go_button_label,
    handle_go_function,
    save_file=True,
    initial_directory=None,
    all_files=False,
  ):
    """GTK4-safe one-entry prompt with a Browse button for open/save paths."""
    window = Gtk.Window()
    window.set_title("Coot")

    vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=8)
    hbox_buttons = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
    entry_row = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
    title_label = Gtk.Label(label=function_label)
    entry = Gtk.Entry()
    browse_button = Gtk.Button(label="Browse...")
    cancel_button = Gtk.Button(label="Cancel")
    go_button = Gtk.Button(label=go_button_label)

    title_label.set_margin_start(12)
    title_label.set_margin_end(12)
    title_label.set_margin_top(12)
    title_label.set_margin_bottom(4)
    title_label.set_xalign(0.0)
    entry_row.set_margin_start(12)
    entry_row.set_margin_end(12)
    entry_row.set_margin_bottom(8)
    entry.set_hexpand(True)
    hbox_buttons.set_margin_start(12)
    hbox_buttons.set_margin_end(12)
    hbox_buttons.set_margin_bottom(12)

    if isinstance(entry_1_default_text, str):
      entry.set_text(entry_1_default_text)

    def close_window(*_args):
      window.destroy()
      return False

    def submit(*_args):
      status = handle_go_function(entry.get_text())
      if status not in (0, False):
        window.destroy()
      return False

    cancel_button.connect("clicked", close_window)
    browse_button.connect(
      "clicked",
      lambda *_args: _populate_entry_from_save_file_chooser(window, entry, entry_1_browse_title)
      if save_file else _populate_entry_from_file_chooser(
        window, entry, entry_1_browse_title, initial_directory, all_files),
    )
    go_button.connect("clicked", submit)
    entry.connect("activate", submit)

    entry_row.append(entry)
    entry_row.append(browse_button)

    vbox.append(title_label)
    vbox.append(entry_row)
    hbox_buttons.append(cancel_button)
    hbox_buttons.append(go_button)
    vbox.append(hbox_buttons)
    window.set_child(vbox)
    window.present()
    entry.grab_focus()

  def generic_title_author_note_entry(
    function_label,
    residue_label,
    title_default_text,
    author_default_text,
    note_default_text,
    go_button_label,
    handle_go_function,
    close_function=None,
  ):
    """Small GTK4 prompt for adding a titled/authored multiline residue note."""
    window = Gtk.Window()
    window.set_title("Coot")
    window.set_default_size(520, 400)

    vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=8)
    hbox_buttons = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
    title_label = Gtk.Label(label=function_label)
    residue_info_label = Gtk.Label(label=residue_label)
    note_title_label = Gtk.Label(label="Title")
    note_title_entry = Gtk.Entry()
    author_label = Gtk.Label(label="Author")
    author_entry = Gtk.Entry()
    note_label = Gtk.Label(label="Note")
    note_scrolled = Gtk.ScrolledWindow()
    note_view = Gtk.TextView()
    note_buffer = note_view.get_buffer()
    cancel_button = Gtk.Button(label="Cancel")
    go_button = Gtk.Button(label=go_button_label)

    title_label.set_margin_start(12)
    title_label.set_margin_end(12)
    title_label.set_margin_top(12)
    title_label.set_margin_bottom(2)
    title_label.set_xalign(0.0)
    residue_info_label.set_margin_start(12)
    residue_info_label.set_margin_end(12)
    residue_info_label.set_margin_bottom(4)
    residue_info_label.set_xalign(0.0)
    note_title_label.set_margin_start(12)
    note_title_label.set_margin_end(12)
    note_title_entry.set_margin_start(12)
    note_title_entry.set_margin_end(12)
    author_label.set_margin_start(12)
    author_label.set_margin_end(12)
    author_entry.set_margin_start(12)
    author_entry.set_margin_end(12)
    note_label.set_margin_start(12)
    note_label.set_margin_end(12)
    note_label.set_margin_top(6)
    note_scrolled.set_margin_start(12)
    note_scrolled.set_margin_end(12)
    note_scrolled.set_margin_bottom(8)
    note_scrolled.set_hexpand(True)
    note_scrolled.set_vexpand(True)
    try:
      note_scrolled.set_has_frame(True)
    except Exception:
      pass
    hbox_buttons.set_margin_start(12)
    hbox_buttons.set_margin_end(12)
    hbox_buttons.set_margin_bottom(12)

    if isinstance(title_default_text, str):
      note_title_entry.set_text(title_default_text)
    if not _annotation_text_or_empty(title_default_text).strip():
      try:
        placeholder_text = _annotation_note_title_placeholder(note_default_text)
        if placeholder_text:
          note_title_entry.set_placeholder_text(placeholder_text)
      except Exception:
        pass
    if isinstance(author_default_text, str):
      author_entry.set_text(author_default_text)
    if isinstance(note_default_text, str):
      note_buffer.set_text(note_default_text)

    note_view.set_wrap_mode(Gtk.WrapMode.WORD_CHAR)
    note_scrolled.set_child(note_view)

    def close_window(*_args):
      if callable(close_function):
        close_function()
      window.destroy()
      return False

    def on_close_request(*_args):
      if callable(close_function):
        close_function()
      return False

    window.connect("close-request", on_close_request)

    def submit(*_args):
      start_iter = note_buffer.get_start_iter()
      end_iter = note_buffer.get_end_iter()
      status = handle_go_function(
        note_title_entry.get_text(),
        author_entry.get_text(),
        note_buffer.get_text(start_iter, end_iter, False),
      )
      if status not in (0, False):
        close_window()
      return False

    cancel_button.connect("clicked", close_window)
    go_button.connect("clicked", submit)
    note_title_entry.connect("activate", lambda *_args: author_entry.grab_focus())
    author_entry.connect("activate", lambda *_args: note_view.grab_focus())

    vbox.append(title_label)
    vbox.append(residue_info_label)
    vbox.append(note_title_label)
    vbox.append(note_title_entry)
    vbox.append(author_label)
    vbox.append(author_entry)
    vbox.append(note_label)
    vbox.append(note_scrolled)
    hbox_buttons.append(cancel_button)
    hbox_buttons.append(go_button)
    vbox.append(hbox_buttons)
    window.set_child(vbox)
    window.present()
    if not author_entry.get_text().strip():
      author_entry.grab_focus()
    elif note_buffer.get_char_count() > 0:
      note_view.grab_focus()
    else:
      note_title_entry.grab_focus()

  def generic_author_note_entry(
    function_label,
    residue_label,
    author_default_text,
    note_default_text,
    go_button_label,
    handle_go_function,
  ):
    return generic_title_author_note_entry(
      function_label,
      residue_label,
      "",
      author_default_text,
      note_default_text,
      go_button_label,
      lambda title_text, author_text, note_text: handle_go_function(author_text, note_text),
    )

  def interesting_things_gui(dialog_name, thing_list):
    """Show a simple GTK4 list of jump targets as labeled buttons."""
    window = Gtk.Window()
    window.set_title("Coot")
    window.set_default_size(520, 420)

    vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=8)
    label = Gtk.Label(label=dialog_name)
    scrolled = Gtk.ScrolledWindow()
    inside_vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6)
    close_button = Gtk.Button(label="Close")

    label.set_margin_start(12)
    label.set_margin_end(12)
    label.set_margin_top(12)
    label.set_margin_bottom(4)
    scrolled.set_margin_start(12)
    scrolled.set_margin_end(12)
    scrolled.set_margin_bottom(8)
    close_button.set_margin_start(12)
    close_button.set_margin_end(12)
    close_button.set_margin_bottom(12)

    # GTK4 will otherwise let the scrolled window collapse to the height of
    # roughly one row, leaving most of the dialog as unused blank space.
    scrolled.set_hexpand(True)
    scrolled.set_vexpand(True)
    inside_vbox.set_valign(Gtk.Align.START)
    scrolled.set_child(inside_vbox)

    def close_window(*_args):
      window.destroy()

    def jump_to_entry(entry):
      _activate_interesting_entry(entry)

    for entry in thing_list:
      if len(entry) < 4:
        continue
      button = Gtk.Button(label=str(entry[0]))
      button.connect(
        "clicked",
        lambda _button, current_entry=entry: jump_to_entry(current_entry),
      )
      inside_vbox.append(button)

    close_button.connect("clicked", close_window)

    vbox.append(label)
    vbox.append(scrolled)
    vbox.append(close_button)
    window.set_child(vbox)
    window.present()

  def navigable_interesting_things_gui(dialog_name, thing_list):
    """Show a flat clickable jump-target list with selection highlighting and Prev/Next."""
    window = Gtk.Window()
    window.set_title("Coot")
    window.set_default_size(520, 420)

    vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=8)
    navigation_row = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
    label = Gtk.Label(label=dialog_name)
    scrolled = Gtk.ScrolledWindow()
    inside_vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6)
    prev_button = Gtk.Button(label="<--Prev")
    next_button = Gtk.Button(label="Next-->")
    close_button = Gtk.Button(label="Close")

    label.set_margin_start(12)
    label.set_margin_end(12)
    label.set_margin_top(12)
    label.set_margin_bottom(4)
    label.set_wrap(True)
    label.set_xalign(0.0)
    scrolled.set_margin_start(12)
    scrolled.set_margin_end(12)
    scrolled.set_margin_bottom(8)
    navigation_row.set_halign(Gtk.Align.CENTER)
    navigation_row.set_margin_start(12)
    navigation_row.set_margin_end(12)
    navigation_row.set_margin_bottom(4)
    close_button.set_margin_start(12)
    close_button.set_margin_end(12)
    close_button.set_margin_bottom(12)
    close_button.set_halign(Gtk.Align.CENTER)

    scrolled.set_hexpand(True)
    scrolled.set_vexpand(True)
    inside_vbox.set_valign(Gtk.Align.START)
    scrolled.set_child(inside_vbox)

    def close_window(*_args):
      window.destroy()

    flat_entry_specs = []
    current_index = {"value": None}
    selected_button = {"value": None}
    category_expanded = {}
    category_widgets = {}

    def set_selected_button(button):
      previous_button = selected_button["value"]
      if previous_button is not None:
        try:
          previous_button.remove_css_class("suggested-action")
        except Exception:
          pass
      selected_button["value"] = button
      if button is not None:
        try:
          button.add_css_class("suggested-action")
        except Exception:
          pass

    def activate_entry(entry_index):
      if not flat_entry_specs:
        return None
      entry_index = entry_index % len(flat_entry_specs)
      current_index["value"] = entry_index
      entry, button = flat_entry_specs[entry_index]
      _activate_interesting_entry(entry)
      set_selected_button(button)
      try:
        button.grab_focus()
      except Exception:
        pass
      return None

    def activate_prev():
      if not flat_entry_specs:
        return None
      if current_index["value"] is None:
        return activate_entry(len(flat_entry_specs) - 1)
      return activate_entry(current_index["value"] - 1)

    def activate_next():
      if not flat_entry_specs:
        return None
      if current_index["value"] is None:
        return activate_entry(0)
      return activate_entry(current_index["value"] + 1)

    for entry in thing_list:
      if len(entry) < 4:
        continue
      button = Gtk.Button(label=str(entry[0]))
      entry_index = len(flat_entry_specs)
      flat_entry_specs.append((entry, button))
      button.connect(
        "clicked",
        lambda _button, idx=entry_index: activate_entry(idx),
      )
      inside_vbox.append(button)

    prev_button.set_sensitive(bool(flat_entry_specs))
    next_button.set_sensitive(bool(flat_entry_specs))
    prev_button.connect("clicked", lambda *_args: activate_prev())
    next_button.connect("clicked", lambda *_args: activate_next())
    close_button.connect("clicked", close_window)

    vbox.append(label)
    vbox.append(scrolled)
    navigation_row.append(prev_button)
    navigation_row.append(next_button)
    vbox.append(navigation_row)
    vbox.append(close_button)
    window.set_child(vbox)
    window.present()

  def categorized_interesting_things_gui(dialog_name, categorized_thing_lists, refresh_function=None, chain_choices=None):
    """Show grouped clickable jump targets in a single scrollable dialog."""
    window = Gtk.Window()
    window.set_title("Coot")
    window.set_default_size(560, 520)
    closed = {"value": False}

    vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=8)
    navigation_row = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
    label = Gtk.Label(label=dialog_name)
    scrolled = Gtk.ScrolledWindow()
    inside_vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=8)
    prev_button = Gtk.Button(label="<--Prev")
    next_button = Gtk.Button(label="Next-->")
    refresh_button = Gtk.Button(label="Refresh list")
    close_button = Gtk.Button(label="Close")

    label.set_margin_start(12)
    label.set_margin_end(12)
    label.set_margin_top(12)
    label.set_margin_bottom(4)
    label.set_wrap(True)
    label.set_xalign(0.0)
    scrolled.set_margin_start(12)
    scrolled.set_margin_end(12)
    scrolled.set_margin_bottom(8)
    navigation_row.set_halign(Gtk.Align.CENTER)
    navigation_row.set_margin_start(12)
    navigation_row.set_margin_end(12)
    navigation_row.set_margin_bottom(4)
    close_button.set_margin_start(12)
    close_button.set_margin_end(12)
    close_button.set_margin_bottom(12)
    close_button.set_halign(Gtk.Align.CENTER)
    refresh_button.set_halign(Gtk.Align.CENTER)

    scrolled.set_hexpand(True)
    scrolled.set_vexpand(True)
    inside_vbox.set_valign(Gtk.Align.START)
    scrolled.set_child(inside_vbox)

    def close_window(*_args):
      closed["value"] = True
      _clear_odd_helical_feature_highlights()
      population_epoch["value"] += 1
      window.destroy()

    def jump_to_entry(entry, category_name=None):
      _activate_interesting_entry(entry, category_name)

    # Keep a flat ordered list for Prev/Next while still rendering grouped rows.
    # `None` means the user has not navigated the list yet.
    flat_entry_specs = []
    current_index = {"value": None}
    selected_button = {"value": None}
    category_expanded = {}
    category_widgets = {}
    population_epoch = {"value": 0}
    glib = _coot_gui_repository_module("GLib")
    filter_state = {"choices": [(None, "All chains")] + list(chain_choices or []),
                    "chain": None, "updating": False, "title": dialog_name,
                    "categories": categorized_thing_lists}
    chain_dropdown = None
    if chain_choices is not None:
      chain_dropdown = Gtk.DropDown.new_from_strings([text for _chain, text in filter_state["choices"]])
      chain_dropdown.set_selected(0)

    def update_chain_choices(choices):
      filter_state["updating"] = True
      try:
        filter_state["choices"] = [(None, "All chains")] + list(choices)
        index = next((i for i, (chain, _text) in enumerate(filter_state["choices"])
                      if chain == filter_state["chain"]), 0)
        filter_state["chain"] = filter_state["choices"][index][0]
        chain_dropdown.set_model(Gtk.StringList.new([text for _chain, text in filter_state["choices"]]))
        chain_dropdown.set_selected(index)
      finally:
        filter_state["updating"] = False

    def set_selected_button(button):
      previous_button = selected_button["value"]
      if previous_button is not None:
        try:
          previous_button.remove_css_class("suggested-action")
        except Exception:
          pass
      selected_button["value"] = button
      if button is not None:
        try:
          button.add_css_class("suggested-action")
        except Exception:
          pass

    def ensure_entry_button(entry_index):
      category_name, entry, button = flat_entry_specs[entry_index]
      if button is None:
        button = Gtk.Button(label=str(entry[0]))
        button.connect("clicked", lambda _button, idx=entry_index: activate_entry(idx))
        flat_entry_specs[entry_index][2] = button
        box = category_widgets[category_name]["box"]
        indices = category_widgets[category_name]["indices"]
        previous = None
        previous_index = entry_index - 1
        while previous_index >= indices[0]:
          previous = flat_entry_specs[previous_index][2]
          if previous is not None:
            break
          previous_index -= 1
        box.insert_child_after(button, previous)
      return button

    def category_toggled(widget, _param, category_name):
      category_expanded[category_name] = widget.get_expanded()
      state = category_widgets[category_name]
      if not widget.get_expanded() or state["pending"]:
        return
      epoch = population_epoch["value"]
      state["pending"] = True
      def populate():
        if epoch != population_epoch["value"] or not widget.get_expanded():
          state["pending"] = False
          return False
        stop = min(state["cursor"] + 40, len(state["indices"]))
        for index in state["indices"][state["cursor"]:stop]:
          ensure_entry_button(index)
        state["cursor"] = stop
        state["pending"] = stop < len(state["indices"])
        return state["pending"]
      if glib is not None:
        glib.idle_add(populate)
      else:
        while populate():
          pass

    def render_categories(current_dialog_name, current_categorized_thing_lists):
      population_epoch["value"] += 1
      if chain_dropdown is not None:
        filter_state["title"] = current_dialog_name
        filter_state["categories"] = current_categorized_thing_lists
        current_categorized_thing_lists = _odd_residue_filter_categories(
          current_categorized_thing_lists, filter_state["chain"])
      label.set_text(current_dialog_name)
      current_index["value"] = None
      set_selected_button(None)
      flat_entry_specs.clear()
      category_widgets.clear()

      child = inside_vbox.get_first_child()
      while child is not None:
        next_child = child.get_next_sibling()
        inside_vbox.remove(child)
        child = next_child

      shown_any_category = False
      for category_name, thing_list in current_categorized_thing_lists:
        if not thing_list:
          continue
        shown_any_category = True

        expander = Gtk.Expander(label=f"{category_name} ({len(thing_list)})")
        expander.connect("notify::expanded", lambda widget, param, name=category_name:
                         category_toggled(widget, param, name))
        category_box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=4)
        expander.set_child(category_box)
        inside_vbox.append(expander)
        category_widgets[category_name] = {"expander": expander, "box": category_box,
                                           "indices": [], "pending": False, "cursor": 0}

        for entry in thing_list:
          if len(entry) < 4:
            continue
          entry_index = len(flat_entry_specs)
          flat_entry_specs.append([category_name, entry, None])
          category_widgets[category_name]["indices"].append(entry_index)
        expander.set_expanded(category_expanded.get(category_name, False))

        separator = Gtk.Separator(orientation=Gtk.Orientation.HORIZONTAL)
        separator.set_margin_top(4)
        separator.set_margin_bottom(2)
        inside_vbox.append(separator)

      if not shown_any_category:
        empty_label = Gtk.Label(label="No results")
        empty_label.set_xalign(0.0)
        inside_vbox.append(empty_label)

      prev_button.set_sensitive(bool(flat_entry_specs))
      next_button.set_sensitive(bool(flat_entry_specs))

    def activate_entry(entry_index):
      if not flat_entry_specs:
        return None
      entry_index = entry_index % len(flat_entry_specs)
      current_index["value"] = entry_index
      category_name, entry, button = flat_entry_specs[entry_index]
      category_widgets[category_name]["expander"].set_expanded(True)
      button = ensure_entry_button(entry_index)
      jump_to_entry(entry, category_name)
      set_selected_button(button)
      try:
        button.grab_focus()
      except Exception:
        pass

    def activate_prev():
      if not flat_entry_specs:
        return None
      if current_index["value"] is None:
        return activate_entry(len(flat_entry_specs) - 1)
      return activate_entry(current_index["value"] - 1)

    def activate_next():
      if not flat_entry_specs:
        return None
      if current_index["value"] is None:
        return activate_entry(0)
      return activate_entry(current_index["value"] + 1)

    def refresh_categories(*_args):
      if not callable(refresh_function):
        return None
      def completed(refreshed_payload):
        if closed["value"]:
          return
        refresh_button.set_sensitive(True)
        if not refreshed_payload:
          return
        refreshed_dialog_name, refreshed_categorized_thing_lists = refreshed_payload[:2]
        if chain_dropdown is not None and len(refreshed_payload) > 2:
          update_chain_choices(refreshed_payload[2])
        render_categories(refreshed_dialog_name, refreshed_categorized_thing_lists)
      start_async = getattr(refresh_function, "start_async", None)
      if callable(start_async):
        refresh_button.set_sensitive(False)
        return start_async(completed)
      return completed(refresh_function())

    def chain_changed(*_args):
      if filter_state["updating"]:
        return
      index = chain_dropdown.get_selected()
      if index < len(filter_state["choices"]):
        filter_state["chain"] = filter_state["choices"][index][0]
        render_categories(filter_state["title"], filter_state["categories"])

    prev_button.connect("clicked", lambda *_args: activate_prev())
    next_button.connect("clicked", lambda *_args: activate_next())
    refresh_button.connect("clicked", refresh_categories)
    close_button.connect("clicked", close_window)
    window.connect("close-request", lambda *_: close_window() or False)
    refresh_button.set_sensitive(callable(refresh_function))

    render_categories(dialog_name, categorized_thing_lists)

    vbox.append(label)
    if chain_dropdown is not None:
      filter_row = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
      filter_row.set_margin_start(12)
      filter_row.set_margin_end(12)
      filter_row.append(Gtk.Label(label="Chain"))
      chain_dropdown.set_hexpand(True)
      filter_row.append(chain_dropdown)
      chain_dropdown.connect("notify::selected", chain_changed)
      vbox.append(filter_row)
    vbox.append(scrolled)
    navigation_row.append(prev_button)
    navigation_row.append(next_button)
    if callable(refresh_function):
      navigation_row.append(refresh_button)
    vbox.append(navigation_row)
    vbox.append(close_button)
    window.set_child(vbox)
    window.present()

  def _interesting_entry_status_text(entry, category_name=None):
    entry_label = str(entry[0]) if entry else ""
    if entry_label and category_name:
      return f"{entry_label} - {_odd_residue_status_suffix(category_name)}"
    return entry_label

  def _activate_interesting_entry(entry, category_name=None):
    """Jump to a chooser entry, using smart residue navigation when available."""
    status_message = _interesting_entry_status_text(entry, category_name)
    if len(entry) >= 5 and isinstance(entry[4], dict):
      navigation = entry[4]
      if navigation.get("type") == "polymer_range":
        _highlight_odd_helical_range(navigation)
        if status_message:
          add_status_bar_text(status_message)
        return _go_to_navigation_residue(
          navigation["mol_id"], navigation["chain_id"], navigation["resno"],
          navigation.get("ins_code", ""), navigation.get("serial_number"))
      _clear_odd_helical_feature_highlights()
      if navigation.get("type") == "polymer_residue":
        if status_message:
          add_status_bar_text(status_message)
        return _go_to_navigation_residue(
          navigation["mol_id"],
          navigation["chain_id"],
          navigation["resno"],
          navigation.get("ins_code", ""),
          navigation.get("serial_number"),
        )
    else:
      _clear_odd_helical_feature_highlights()
    if status_message:
      add_status_bar_text(status_message)
    set_rotation_centre(entry[1], entry[2], entry[3])
    return False

  def action_button_dialog(dialog_name, button_list, close_on_click=True):
    """GTK4-safe scrolling list of action buttons."""
    window = Gtk.Window()
    window.set_title("Coot")
    window.set_default_size(520, 420)

    vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=8)
    label = Gtk.Label(label=dialog_name)
    scrolled = Gtk.ScrolledWindow()
    inside_vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6)
    close_button = Gtk.Button(label="Close")

    label.set_margin_start(12)
    label.set_margin_end(12)
    label.set_margin_top(12)
    label.set_margin_bottom(4)
    scrolled.set_margin_start(12)
    scrolled.set_margin_end(12)
    scrolled.set_margin_bottom(8)
    close_button.set_margin_start(12)
    close_button.set_margin_end(12)
    close_button.set_margin_bottom(12)

    scrolled.set_hexpand(True)
    scrolled.set_vexpand(True)
    inside_vbox.set_valign(Gtk.Align.START)
    scrolled.set_child(inside_vbox)

    def close_window(*_args):
      window.destroy()

    for button_label, callback in button_list:
      button = Gtk.Button(label=str(button_label))

      def on_click(_button, action=callback):
        action()
        if close_on_click:
          window.destroy()

      button.connect("clicked", on_click)
      inside_vbox.append(button)

    close_button.connect("clicked", close_window)

    vbox.append(label)
    vbox.append(scrolled)
    vbox.append(close_button)
    window.set_child(vbox)
    window.present()

else:

  def coot_toolbar_button(*_args, **_kwargs):
    return None

  def attach_module_menu_button(_module_name):
    return None

  def generic_button_dialog(dialog_name, _button_list):
    info_dialog(_gui_unavailable_message(dialog_name))

  def generic_multiple_entries_with_check_button(
    _entry_info_list, _check_button_info, go_button_label, _handle_go_function
  ):
    info_dialog(_gui_unavailable_message(go_button_label))

  def generic_single_entry(
    function_label, _entry_1_default_text, _go_button_label, _handle_go_function
  ):
    info_dialog(_gui_unavailable_message(function_label))

  def generic_confirm_dialog(
    function_label,
    _message_text,
    _cancel_button_label,
    _handle_cancel_function,
    _confirm_button_label,
    _handle_confirm_function,
    handle_close_function=None,
  ):
    if handle_close_function is not None:
      handle_close_function()
    info_dialog(_gui_unavailable_message(function_label))

  def generic_double_entry(
    function_label,
    _entry_1_label,
    _entry_1_default_text,
    _entry_2_label,
    _entry_2_default_text,
    _go_button_label,
    _handle_go_function,
  ):
    info_dialog(_gui_unavailable_message(function_label))

  def generic_double_entry_with_file_browse(
    function_label,
    _entry_1_label,
    _entry_1_default_text,
    _entry_2_label,
    _entry_2_default_text,
    _entry_2_browse_title,
    _go_button_label,
    _handle_go_function,
  ):
    info_dialog(_gui_unavailable_message(function_label))

  def generic_single_entry_with_file_browse(
    function_label,
    _entry_1_default_text,
    _entry_1_browse_title,
    _go_button_label,
    _handle_go_function,
    save_file=True,
    initial_directory=None,
    all_files=False,
  ):
    info_dialog(_gui_unavailable_message(function_label))

  def generic_title_author_note_entry(
    function_label,
    _residue_label,
    _title_default_text,
    _author_default_text,
    _note_default_text,
    _go_button_label,
    _handle_go_function,
    close_function=None,
  ):
    info_dialog(_gui_unavailable_message(function_label))
    if callable(close_function):
      close_function()

  def generic_author_note_entry(
    function_label,
    _residue_label,
    _author_default_text,
    _note_default_text,
    _go_button_label,
    _handle_go_function,
  ):
    info_dialog(_gui_unavailable_message(function_label))

  def interesting_things_gui(dialog_name, thing_list):
    labels = [str(entry[0]) for entry in thing_list if len(entry) >= 4]
    message = dialog_name
    if labels:
      message += "\n\n" + "\n".join(labels)
    info_dialog(message)

  def navigable_interesting_things_gui(dialog_name, thing_list):
    interesting_things_gui(dialog_name, thing_list)

  def categorized_interesting_things_gui(dialog_name, categorized_thing_lists, refresh_function=None, chain_choices=None):
    labels = []
    for category_name, thing_list in categorized_thing_lists:
      if not thing_list:
        continue
      labels.append(f"{category_name} ({len(thing_list)})")
      labels.extend(str(entry[0]) for entry in thing_list if len(entry) >= 4)
    message = dialog_name
    if labels:
      message += "\n\n" + "\n".join(labels)
    info_dialog(message)

  def action_button_dialog(dialog_name, button_list, _close_on_click=True):
    labels = [str(label) for label, _callback in button_list]
    message = dialog_name
    if labels:
      message += "\n\n" + "\n".join(labels)
    info_dialog(message)


  def add_simple_coot_menu_menuitem(_menu, _menu_item_label, _activate_function):
    return None

  def add_module_cryo_em_gui():
    return None

  def add_module_refine():
    return None


def add_key_binding(name, key, thunk):
  if not hasattr(coot, "add_key_binding_gtk4_py"):
    _show_startup_diagnostic_once(
      "missing_gtk4_keybinding_api",
      "This Coot build does not expose add_key_binding_gtk4_py, so coot_trimmings "
      "custom keybindings are disabled. Menus and non-keybinding helpers may still load.",
    )
    return None
  ctrl_key = 0
  key_value = key
  if isinstance(key, str) and key.startswith("Control_"):
    ctrl_key = 1
    key_value = key[len("Control_"):]
  if isinstance(key_value, str) and len(key_value) > 1:
    named_keys = {"Tab": 65289, "Return": 65293, "Escape": 65307,
                  "Left": 65361, "Up": 65362, "Right": 65363, "Down": 65364,
                  "BackSpace": 65288, "Delete": 65535, "space": 32}
    if key_value not in named_keys:
      raise ValueError(f"Unknown named key: {key_value}")
    key_value = named_keys[key_value]

  def wrapped_thunk():
    try:
      return thunk()
    except Exception:
      print("coot_trimmings keybinding failure:", name, "key=", key_value)
      traceback.print_exc()
      return None

  REGISTERED_KEYBINDING_CALLBACKS.append(wrapped_thunk)
  coot.add_key_binding_gtk4_py(key_value, ctrl_key, wrapped_thunk, name)
  key_id = ord(key_value) if isinstance(key_value, str) else key_value
  REGISTERED_KEYBINDINGS[(ctrl_key, key_id)] = (key, name)
  globals().setdefault("REGISTERED_KEYBINDING_ACTIONS", {})[(ctrl_key, key_id)] = {
      "key": key, "name": name, "callback": wrapped_thunk}


def _custom_key_binding_label(key):
  if isinstance(key, str):
    return "Ctrl+" + key[len("Control_"):] if key.startswith("Control_") else key
  named = {65289: "Tab", 65293: "Return", 65307: "Escape", 65361: "Left",
           65362: "Up", 65363: "Right", 65364: "Down", 65288: "BackSpace",
           65535: "Delete", 32: "Space", 65421: "Keypad Enter", 65451: "Keypad +", 65453: "Keypad -"}
  if 65470 <= key <= 65504:
    return f"F{key - 65469}"
  if 65456 <= key <= 65465:
    return f"Keypad {key - 65456}"
  return named.get(key, chr(key) if 33 <= key <= 126 else f"Key {key}")


def _parse_native_key_bindings(text, keyval_from_name):
  bindings = {}
  for line in text.splitlines():
    match = re.match(r"^binding:\s*(Ctrl\s+)?(\S+)\s+->\s*(.*?)\s+BUILT_IN\s*$", line)
    if not match:
      continue
    ctrl, symbol, name = match.groups()
    if symbol.startswith("gdk-symbol-") and symbol[11:].isdigit():
      key_id = int(symbol[11:])
    elif len(symbol) == 1:
      key_id = ord(symbol)
    else:
      key_id = keyval_from_name({"KP_plus": "KP_Add", "KP_minus": "KP_Subtract"}.get(symbol, symbol))
    if not key_id or key_id == 0xffffff or not name:
      continue
    identity = (int(bool(ctrl)), key_id)
    if identity in REGISTERED_KEYBINDINGS:
      continue
    key = _custom_key_binding_label(key_id)
    bindings[identity] = {"key": "Ctrl+" + key if ctrl else key, "name": name, "source": "Coot"}
  return list(bindings.values())


def _read_native_key_bindings():
  """Read the live native table only when the reference panel is opened."""
  printer = getattr(coot, "print_key_bindings", None)
  if not callable(printer):
    return [], "Native shortcuts are unavailable in this Coot build."
  try:
    gdk = _coot_gui_repository_module("Gdk")
    # The API writes C++ std::cout, so redirect_stdout cannot capture it.
    # Its lines are flushed with std::endl; restore fd 1 even on API failure.
    sys.stdout.flush()
    with tempfile.TemporaryFile() as output:
      previous = os.dup(1)
      try:
        os.dup2(output.fileno(), 1)
        printer()
      finally:
        os.dup2(previous, 1)
        os.close(previous)
      output.seek(0)
      text = output.read().decode("utf-8", errors="replace")
    for line in text.splitlines():
      if not line.startswith("binding:"):
        print(line)
    return _parse_native_key_bindings(text, gdk.keyval_from_name), None
  except Exception as error:
    print(f"Could not read native key bindings: {error}")
    return [], "Native shortcuts could not be read; showing Custom shortcuts only."


def _key_binding_reference_groups(native, query=""):
  groups = {name: [] for name in ("Build / Edit", "Refinement", "Maps", "Navigation / Display", "Other")}
  items = [{"key": _custom_key_binding_label(key), "name": name, "source": "Custom"}
           for key, name in REGISTERED_KEYBINDINGS.values()] + native
  terms = query.strip().casefold().split()
  for item in items:
    key, action, source = item["key"], item["name"], item["source"]
    if not all(source.casefold() == term if term in ("custom", "coot") else
               term in f"{key} {action}".casefold() for term in terms):
      continue
    group = _key_binding_category(action)
    groups[group].append((key, action, source))
  return [(name, rows) for name, rows in groups.items() if rows]


def _key_binding_category(action):
  lower = action.casefold()
  if any(word in lower for word in ("map", "density", "contour", "sigma", "surface")):
    return "Maps"
  if any(word in lower for word in ("representation", "display", "symmetry", "clipping", "ncs",
                                   "distances", "measure", "pending picks", "residue in chain",
                                   "zoom", "move forward", "move backward", "label", "focus command")):
    return "Navigation / Display"
  if any(word in lower for word in ("refine", "fit", "restraint")):
    return "Refinement"
  if any(word in lower for word in ("place", "copy", "paste", "delete", "sidechain", "water", "mutate",
                                   "terminal", "terminus", "peptide", "rotamer", "undo", "redo", "save")):
    return "Build / Edit"
  return "Other"


def _custom_key_binding_groups(query=""):
  groups = {name: [] for name in ("Build / Edit", "Refinement", "Maps", "Navigation / Display", "Other")}
  query = query.strip().casefold()
  for key, action in REGISTERED_KEYBINDINGS.values():
    key = _custom_key_binding_label(key)
    if query and query not in f"{key} {action}".casefold():
      continue
    group = _key_binding_category(action)
    groups[group].append((key, action))
  return [(name, rows) for name, rows in groups.items() if rows]


def custom_key_bindings_dialog():
  """A read-only reference to Custom shortcuts and remaining native bindings."""
  global _CUSTOM_KEY_BINDINGS_WINDOW
  if Gtk is None:
    info_dialog(_gui_unavailable_message("Custom Key Bindings"))
    return None
  if _CUSTOM_KEY_BINDINGS_WINDOW is not None and _CUSTOM_KEY_BINDINGS_WINDOW.get_visible():
    _CUSTOM_KEY_BINDINGS_WINDOW.key_bindings_widgets["reload"]()
    _CUSTOM_KEY_BINDINGS_WINDOW.present()
    return _CUSTOM_KEY_BINDINGS_WINDOW
  window = Gtk.Window(title="Key Bindings: Custom and Coot")
  window.set_default_size(600, 500)
  box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=8)
  for edge in ("start", "end", "top", "bottom"):
    getattr(box, "set_margin_" + edge)(12)
  hint = Gtk.Label(label="Keys are case-sensitive; uppercase letters require Shift.\nUse shortcuts with focus in the model view. Overridden native keys are omitted.")
  hint.set_xalign(0)
  hint.set_wrap(True)
  box.append(hint)
  availability = Gtk.Label()
  availability.set_xalign(0)
  availability.set_wrap(True)
  availability.add_css_class("dim-label")
  box.append(availability)
  search = Gtk.SearchEntry()
  search.set_placeholder_text("Find a key or action; filter with custom / coot")
  box.append(search)
  scroll = Gtk.ScrolledWindow()
  scroll.set_policy(Gtk.PolicyType.NEVER, Gtk.PolicyType.AUTOMATIC)
  scroll.set_hexpand(True)
  scroll.set_vexpand(True)
  content = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=12)
  scroll.set_child(content)
  box.append(scroll)
  footer = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
  count = Gtk.Label()
  count.set_xalign(0)
  count.set_hexpand(True)
  footer.append(count)
  close = Gtk.Button(label="Close")
  footer.append(close)
  box.append(footer)
  rows = []
  entries, native = [], []
  def refresh(*_args):
    child = content.get_first_child()
    while child is not None:
      content.remove(child)
      child = content.get_first_child()
    rows.clear()
    entries.clear()
    groups = _key_binding_reference_groups(native, search.get_text())
    for title, bindings in groups:
      section = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=4)
      heading = Gtk.Label(label=title)
      heading.set_xalign(0)
      heading.add_css_class("heading")
      section.append(heading)
      grid = Gtk.Grid(column_spacing=16, row_spacing=4)
      for index, (key, action, source) in enumerate(bindings):
        key_label = Gtk.Label(label=key)
        key_label.add_css_class("monospace")
        key_label.set_xalign(0)
        key_label.set_width_chars(10)
        key_label.set_selectable(True)
        action_label = Gtk.Label(label=action)
        action_label.set_xalign(0)
        action_label.set_wrap(True)
        action_label.set_hexpand(True)
        action_label.set_selectable(True)
        grid.attach(key_label, 0, index, 1, 1)
        grid.attach(action_label, 1, index, 1, 1)
        source_label = Gtk.Label(label=source)
        source_label.add_css_class("dim-label")
        source_label.set_xalign(0)
        grid.attach(source_label, 2, index, 1, 1)
        rows.append((key, action))
        entries.append((key, action, source))
      section.append(grid)
      content.append(section)
    total = len(REGISTERED_KEYBINDINGS) + len(native)
    count.set_text(f"{len(rows)} / {total} shortcuts" if search.get_text().strip() else f"{total} shortcuts")
    if not rows:
      content.append(Gtk.Label(label="No matching shortcuts." if total else "No custom shortcuts registered in this build."))
  def reload():
    bindings, note = _read_native_key_bindings()
    native[:] = bindings
    availability.set_text(note or "")
    availability.set_visible(bool(note))
    refresh()
  def close_window(*_args):
    global _CUSTOM_KEY_BINDINGS_WINDOW
    _CUSTOM_KEY_BINDINGS_WINDOW = None
    window.destroy()
    return True
  search.connect("search-changed", refresh)
  close.connect("clicked", close_window)
  window.connect("close-request", close_window)
  window.set_child(box)
  window.key_bindings_widgets = {"search": search, "content": content, "count": count,
                                "close": close, "rows": rows, "entries": entries, "native": native,
                                "availability": availability, "refresh": refresh, "reload": reload}
  _CUSTOM_KEY_BINDINGS_WINDOW = window
  reload()
  window.present()
  return window

if not CARBOHYDRATE_GUI_AVAILABLE:
  def add_module_carbohydrate_gui():
    return None

_apply_startup_settings()

# ============================================================================
# Keybindings (user-editable)
# ============================================================================
#
# This is the main place to tweak keyboard shortcuts. The bindings live here
# rather than at the very top so they can refer directly to the helpers they
# call.
#

#place helix with prosmart alpha helix restraints and depict in rainbow
add_key_binding("Place helix here","h",
lambda: place_helix_with_restraints())

#Measure distance
add_key_binding("Measure distance","m",
lambda: start_measure_distance())

#Go to nearest map peak around the rotation centre
add_key_binding("Go to nearest density peak","b",
lambda: go_to_nearest_density_peak())

#Resample active map to 0.5 A/pixel and restyle
add_key_binding("Resample active map to 0.5 A/pixel","C",
lambda: resample_active_map_for_em_half_angstrom())

#Smart-copy active ligand/ion/water
add_key_binding("Smart copy active non-polymer residue","c",
lambda: smart_copy_active_non_polymer_residue())

#Paste last smart-copied ligand/ion/water at pointer
add_key_binding("Smart paste copied non-polymer residue","v",
lambda: smart_paste_copied_non_polymer_residue())

#Toggle global display of map
add_key_binding("Toggle global view of map","G",
lambda: toggle_global_map_view())

#Darken/Brighten current scrollable map
add_key_binding("Darken current map",",",
lambda: darken_scrollable_map())
add_key_binding("Brighten current map",".",
lambda: brighten_scrollable_map())
add_key_binding("Set map contour in sigma","L",
lambda: set_map_level_quickly())

#Generate smart extra distance restraints on active model
add_key_binding("Generate smart local extra restraints","g",
lambda: generate_smart_local_extra_restraints())


#Quicksave active mol (overwrite orig)
add_key_binding("Save and overwrite active model","Q",
lambda: quicksave_active())

#Refine zone (click two atoms)
add_key_binding("Refine zone","A",
lambda: refine_residues_click())

#Flip peptide.
add_key_binding("Flip peptide","q",
lambda: flip_active_peptide())

#Local cylinder refinement around active residue
add_key_binding("Local cylinder refine","a",
lambda: auto_refine())

#Jiggle fit active res
add_key_binding("Jiggle Fit","J",
lambda: jiggle_fit_active_non_polymer_residue())

#Clear pending picks
add_key_binding("Clear Pending Picks","Tab",
lambda: clear_pending_custom_picks())

#Toggle environment distances
add_key_binding("Toggle environment distances","D",
lambda: toggle_env_dist()) 

#Delete active residue
add_key_binding("Delete this residue","X",
lambda: using_active_atom(delete_residue,"aa_imol","aa_chain_id","aa_res_no","aa_ins_code"))

#Delete sidechain
add_key_binding("Kill Sidechain","K",
lambda: using_active_atom(delete_residue_sidechain,"aa_imol","aa_chain_id","aa_res_no","aa_ins_code",0))

#Fill sidechain
add_key_binding("Fill Sidechain","k",
lambda: using_active_atom(fill_partial_residue,"aa_imol","aa_chain_id","aa_res_no","aa_ins_code"))

#place water without refinement
add_key_binding("Add Water","w",
lambda: place_water_in_active_molecule())

#place water with refinement
add_key_binding("Add Water +","W",
lambda: add_water_and_refine())

#add terminal residue
add_key_binding("Add terminal residue","y",
lambda: add_term_shortcut())
 
#add terminal residue
add_key_binding("Cycle terminus phi","Y",
lambda: cycle_residue_phi())

#add terminal residue
add_key_binding("Cycle terminus psi","T",
lambda: cycle_residue_psi())

 
#Refine active residue
add_key_binding("Refine Triple","r",
lambda: key_binding_refine_triple())

#Undo Symm view
add_key_binding("Undo Symmetry View", "V",
lambda: undo_symmetry_view_safe())

#Cycle through rotamers for active reidue with 'R"
add_key_binding("Cycle rotamers","R",
lambda: cycle_rotamers())

#Undo function for keybinding. Undoes last change to active mol.
add_key_binding("Undo","z",
lambda: undo_visible())

#Redo function for keybinding. redoes last change to active mol.
add_key_binding("Redo","x",
lambda: redo_visible())


add_key_binding("Cycle representation mode forward","[",
lambda: cycle_rep_up_current())

add_key_binding("Cycle representation mode back","]",
lambda: cycle_rep_down_current())

add_key_binding("Cycle  symm representation mode forward","{",
lambda: cycle_symm_up_current())

add_key_binding("Cycle  symm representation mode back","}",
lambda: cycle_symm_down_current())

add_key_binding("Toggle map display","`",
lambda: toggle_map_display())

add_key_binding("Toggle mol display","/",
lambda: toggle_mol_display())

#Clear distances/labels
add_key_binding("Clear distances and labels","Z",
lambda: clear_distances_and_labels())

#Bind next_res() and prev_res() to ">" and "<"
add_key_binding("Next residue in chain",">",
lambda: next_res())
add_key_binding("Prev residue in chain","<",
lambda: prev_res())

#The nine key bindings elow allow easy setting of map
#level by rmsd - shift + any single digit integer
#sets the currently scrollable map to that level
# in rmsd. Useful when on a laptop with touchpad,
#when changing the contour using the scrollwheel is
#not practical and using +/- is too slow.
add_key_binding("Map to 1 sigma","!",
lambda: set_current_map_sigma(1))

add_key_binding("Map to 2 sigma","@",
lambda: set_current_map_sigma(2))

add_key_binding("Map to 3 sigma","#",
lambda: set_current_map_sigma(3))

add_key_binding("Map to 4 sigma","$",
lambda: set_current_map_sigma(4))

add_key_binding("Map to 5 sigma","%",
lambda: set_current_map_sigma(5))

add_key_binding("Map to 6 sigma","^",
lambda: set_current_map_sigma(6))

add_key_binding("Map to 7 sigma","&",
lambda: set_current_map_sigma(7))

add_key_binding("Map to 8 sigma","*",
lambda: set_current_map_sigma(8))

add_key_binding("Map to 9 sigma","(",
lambda: set_current_map_sigma(9))

add_key_binding("Optimize local contour by thresholded Q-residue",")",
lambda: optimize_local_threshold_by_thresholded_q_current_view())

add_key_binding("Map plus 0.1 sigma","+",
lambda: step_current_map_coarse_up())

add_key_binding("Map minus 0.1 sigma","_",
lambda: step_current_map_coarse_down())

add_key_binding("Increase map radius","'",
lambda: increase_map_radius())

add_key_binding("Decrease map radius",";",
lambda: decrease_map_radius())

add_key_binding("Increase active map surface opacity","\"",
lambda: increase_active_map_surface_opacity())

add_key_binding("Decrease active map surface opacity",":",
lambda: decrease_active_map_surface_opacity())

#Undisplay all models except the active one.
#If only one model is displayed, cycle through
#all available models.

add_key_binding("Display only the active model","?",
lambda: display_only_active())  



#Undisplay all maps except the active one.
#If only one map is displayed, cycle through
#all available models.
add_key_binding("Display only the active map","~",
lambda: display_only_active_map())

#Go to equivalent residue on ncs master chain

def goto_ncs_master():
  residue = _active_residue_or_status()
  if not residue:
    return None
  mol_id=residue[0]
  resno=residue[2]
  atom_name=residue[4]
  ncs_ch_id=ncs_master_chain_id(mol_id)
  if not ncs_ch_id:
    add_status_bar_text("No NCS master chain found")
    return None
  set_go_to_atom_molecule(mol_id)
  set_go_to_atom_chain_residue_atom_name(ncs_ch_id,resno,atom_name)
add_key_binding("Go to NCS master chain","O",
lambda: goto_ncs_master())

add_key_binding("Narrow clipping slab","-",
lambda: narrow_clipping_symmetric())

add_key_binding("Widen clipping slab","=",
lambda: widen_clipping_symmetric())

  
#****Misc. functions (for keybindings and scripting****
def toggle_high_contrast_mode():
  """Toggle a simple ambient-only lighting preset for all visible models."""
  global MODEL_AMBIENT_LIGHTING_ENABLED
  global MODEL_PRE_HIGH_CONTRAST_GL_LIGHTING_STATE
  global MODEL_PRE_HIGH_CONTRAST_BOND_THICKNESS
  global MODEL_HIGH_CONTRAST_MOLECULES
  current_model_molecules = set(model_molecule_list())
  needs_reapply = (
    MODEL_AMBIENT_LIGHTING_ENABLED
    and bool(current_model_molecules)
    and not current_model_molecules.issubset(MODEL_HIGH_CONTRAST_MOLECULES)
  )
  if MODEL_AMBIENT_LIGHTING_ENABLED and not needs_reapply:
    ambient = MODEL_NORMAL_AMBIENT
    diffuse = MODEL_NORMAL_DIFFUSE
    specular = MODEL_NORMAL_SPECULAR
    use_variable_bonds = STARTUP_USE_VARIABLE_BOND_THICKNESS
    default_bond_thickness = MODEL_PRE_HIGH_CONTRAST_BOND_THICKNESS
    if default_bond_thickness is None:
      default_bond_thickness = STARTUP_DEFAULT_BOND_THICKNESS
    set_do_GL_lighting(1 if MODEL_PRE_HIGH_CONTRAST_GL_LIGHTING_STATE else 0)
    status_message = "Set all models to normal lighting"
    MODEL_HIGH_CONTRAST_MOLECULES = set()
  else:
    if not MODEL_AMBIENT_LIGHTING_ENABLED:
      MODEL_PRE_HIGH_CONTRAST_GL_LIGHTING_STATE = do_GL_lighting_state()
      MODEL_PRE_HIGH_CONTRAST_BOND_THICKNESS = get_default_bond_thickness()
    set_do_GL_lighting(1)
    ambient = MODEL_AMBIENT_ONLY_AMBIENT
    diffuse = MODEL_AMBIENT_ONLY_DIFFUSE
    specular = MODEL_AMBIENT_ONLY_SPECULAR
    use_variable_bonds = 0
    default_bond_thickness = HIGH_CONTRAST_BOND_THICKNESS
    status_message = (
      "Reapplied high-contrast ambient lighting to current models"
      if needs_reapply
      else "Set all models to high-contrast ambient lighting"
    )
  set_use_variable_bond_thickness(use_variable_bonds)
  set_default_bond_thickness(default_bond_thickness)
  for imol in current_model_molecules:
    set_model_material_ambient(imol, ambient[0], ambient[1], ambient[2], ambient[3])
    set_model_material_diffuse(imol, diffuse[0], diffuse[1], diffuse[2], diffuse[3])
    set_model_material_specular(imol, specular[0], specular[1])
    set_bond_thickness(imol, default_bond_thickness)
  if needs_reapply:
    MODEL_HIGH_CONTRAST_MOLECULES = set(current_model_molecules)
  else:
    MODEL_AMBIENT_LIGHTING_ENABLED = not MODEL_AMBIENT_LIGHTING_ENABLED
    if MODEL_AMBIENT_LIGHTING_ENABLED:
      MODEL_HIGH_CONTRAST_MOLECULES = set(current_model_molecules)
  add_status_bar_text(status_message)


def _display_preferences_path():
  """Return the companion preferences file beside this startup script."""
  return os.path.join(_coot_trimmings_dir(), DISPLAY_PREFERENCES_FILENAME)


def _load_display_defaults():
  """Load display-default booleans from the companion Python preferences file."""
  defaults = {"high_contrast_models": False, "ambient_maps": False, "black_background": False,
              "refinement_ramachandran_overlay": True, "refinement_validation_spheres": None,
              "refinement_contact_dots": None}
  path = _display_preferences_path()
  namespace = {}
  try:
    with open(path, "r", encoding="utf-8") as handle:
      source = handle.read()
    exec(compile(source, path, "exec"), {"__builtins__": {}}, namespace)
    defaults["high_contrast_models"] = bool(namespace.get("DEFAULT_HIGH_CONTRAST_MODELS", False))
    defaults["ambient_maps"] = bool(namespace.get("DEFAULT_AMBIENT_MAPS", False))
    defaults["black_background"] = bool(namespace.get("DEFAULT_BLACK_BACKGROUND", False))
    defaults["refinement_ramachandran_overlay"] = bool(namespace.get("DEFAULT_REFINEMENT_RAMACHANDRAN_OVERLAY", True))
    for name in ("refinement_validation_spheres", "refinement_contact_dots"):
      value = namespace.get("DEFAULT_" + name.upper())
      if value is not None:
        defaults[name] = bool(value)
  except FileNotFoundError:
    pass
  except Exception as error:
    print(f"WARNING:: could not read display preferences from {path}: {error}")
  return defaults


def _save_display_default(name, enabled=True):
  """Persist display defaults as simple Python beside coot_trimmings.py."""
  defaults = _load_display_defaults()
  defaults[name] = bool(enabled)
  path = _display_preferences_path()
  source = (
    '# Auto-generated by coot_trimmings.py. You may edit these values manually.\n'
    f'DEFAULT_HIGH_CONTRAST_MODELS = {defaults["high_contrast_models"]!r}\n'
    f'DEFAULT_AMBIENT_MAPS = {defaults["ambient_maps"]!r}\n'
    f'DEFAULT_BLACK_BACKGROUND = {defaults["black_background"]!r}\n'
    f'DEFAULT_REFINEMENT_RAMACHANDRAN_OVERLAY = {defaults["refinement_ramachandran_overlay"]!r}\n'
  )
  for key in ("refinement_validation_spheres", "refinement_contact_dots"):
    if defaults[key] is not None:
      source += f'DEFAULT_{key.upper()} = {defaults[key]!r}\n'
  try:
    # Atomic replacement avoids leaving a partial preferences file if writing fails.
    directory = os.path.dirname(path) or "."
    fd, temporary_path = tempfile.mkstemp(prefix=".coot_trimmings_preferences_", suffix=".py", dir=directory)
    try:
      with os.fdopen(fd, "w", encoding="utf-8") as handle:
        handle.write(source)
      os.replace(temporary_path, path)
    except Exception:
      try:
        os.unlink(temporary_path)
      except Exception:
        pass
      raise
  except Exception as error:
    info_dialog(
      "Could not save the display default for future sessions.\n\n"
      f"Preferences file: {path}\n\n{error}"
    )
    return False
  return True


def reset_map_model_defaults():
  """Remove persisted display defaults; Coot defaults return after restart."""
  path = _display_preferences_path()
  try:
    os.unlink(path)
  except FileNotFoundError:
    add_status_bar_text("Display defaults are already reset; restart Coot to use Coot defaults")
    return True
  except Exception as error:
    info_dialog(
      "Could not reset the persisted display defaults.\n\n"
      f"Preferences file: {path}\n\n{error}"
    )
    return False
  add_status_bar_text("Display defaults reset; restart Coot to restore Coot defaults")
  return True


def _set_black_background(show_dialog=True):
  try:
    coot.set_background_colour(0.0, 0.0, 0.0)
  except Exception as error:
    if show_dialog:
      info_dialog(f"Could not set a black background in this Coot build.\n\n{error}")
    else:
      print(f"WARNING:: could not restore black background: {error}")
    return False
  return True


def default_background_to_black():
  """Apply black now and restore it during subsequent startups."""
  if not _set_black_background() or not _save_display_default("black_background", True):
    return False
  add_status_bar_text("Black background enabled for this and future Coot sessions")
  return True


def _set_refinement_ramachandran_overlay(enabled, show_dialog=True):
  """Apply and track the plot visibility; silently skip unsupported startup APIs."""
  global REFINEMENT_RAMACHANDRAN_OVERLAY_ENABLED
  setter = getattr(coot, "set_draw_gl_ramachandran_plot_during_refinement", None)
  if not callable(setter):
    if show_dialog:
      info_dialog("This Coot build cannot toggle the Ramachandran overlay during refinement.")
    return False
  try:
    setter(int(bool(enabled)))
  except Exception as error:
    if show_dialog:
      info_dialog(f"Could not change the refinement Ramachandran overlay.\n\n{error}")
    else:
      print(f"WARNING:: could not restore refinement Ramachandran overlay: {error}")
    return False
  REFINEMENT_RAMACHANDRAN_OVERLAY_ENABLED = bool(enabled)
  return True


def _refinement_validation_display_controls(name):
  """Return native getters/setters for the requested refinement-only display."""
  if name == "refinement_validation_spheres":
    label = "Refinement Rama/rotamer spheres"
    functions = [("get_draw_moving_atoms_rama_markup_state", "set_draw_moving_atoms_rama_markup"),
                 ("get_draw_moving_atoms_rota_markup_state", "set_draw_moving_atoms_rota_markup")]
  elif name == "refinement_contact_dots":
    label = "Refinement contact dots"
    functions = [("get_do_coot_probe_dots_during_refine", "set_do_coot_probe_dots_during_refine")]
  else:
    raise ValueError("Unknown refinement display setting: " + name)
  controls = [(getattr(coot, getter, None), getattr(coot, setter, None)) for getter, setter in functions]
  if not all(callable(function) for pair in controls for function in pair):
    raise RuntimeError(label + " cannot be toggled in this Coot build.")
  return label, controls


def _set_refinement_validation_display(name, enabled, show_dialog=True):
  """Restore saved display flags without changing refinement restraints."""
  previous = []
  try:
    _label, controls = _refinement_validation_display_controls(name)
    previous = [(setter, int(bool(getter()))) for getter, setter in controls]
    for _getter, setter in controls:
      setter(int(bool(enabled)))
  except Exception as error:
    for setter, value in previous:
      try:
        setter(value)
      except Exception as restore_error:
        print(f"WARNING:: could not restore refinement display flag: {restore_error}")
    if show_dialog:
      info_dialog(f"Could not change refinement display.\n\n{error}")
    return False
  return True


def _toggle_refinement_validation_display(name):
  try:
    label, controls = _refinement_validation_display_controls(name)
    # A mixed Rama/rotamer state switches both off; the next click enables both.
    enabled = not any(bool(getter()) for getter, _setter in controls)
  except Exception as error:
    info_dialog(f"Could not change refinement display.\n\n{error}")
    return None
  if not _set_refinement_validation_display(name, enabled):
    return None
  saved = _save_display_default(name, enabled)
  add_status_bar_text(label + ": " + ("on" if enabled else "off") +
                      (" (saved for future sessions)" if saved else " (this session only; could not save)"))
  return enabled


def toggle_refinement_validation_spheres():
  return _toggle_refinement_validation_display("refinement_validation_spheres")


def toggle_refinement_contact_dots():
  # Temporarily disabled pending Coot's contact-dot display fix. Keep saved choices.
  add_status_bar_text("Refinement contact-dot toggle is temporarily disabled pending a Coot fix")
  return None


def _trackpad_preferences_path():
  return os.path.join(_coot_trimmings_dir(), TRACKPAD_PREFERENCES_FILENAME)


def _trackpad_normalize_preferences(data):
  if not isinstance(data, dict) or data.get("version", 1) != 1:
    raise ValueError("Unsupported trackpad preferences")
  result = {"version": 1, "startup": False, "reverse": True, "sensitivity": {}}
  for name in ("startup", "reverse"):
    value = data.get(name, result[name])
    if type(value) is not bool:
      raise ValueError(f"Invalid trackpad {name} setting")
    result[name] = value
  values = data.get("sensitivity", {})
  if not isinstance(values, dict):
    raise ValueError("Invalid trackpad sensitivities")
  for name in ("zoom", "rotate", "translate", "contour", "clipping"):
    value = values.get(name, 2.0)
    if type(value) not in (float, int) or not math.isfinite(value) or not 0.1 <= value <= 3.0:
      raise ValueError(f"Invalid trackpad {name} sensitivity")
    result["sensitivity"][name] = float(value)
  return result


def _save_trackpad_preferences(startup=None):
  state = _TRACKPAD_STATE
  temporary = None
  try:
    data = _trackpad_normalize_preferences({
      "startup": state["startup"] if startup is None else startup,
      "reverse": state["reverse"], "sensitivity": state["sensitivity"]})
    path = _trackpad_preferences_path()
    fd, temporary = tempfile.mkstemp(prefix=".coot_trackpad_", suffix=".json", dir=os.path.dirname(path))
    with os.fdopen(fd, "w", encoding="utf-8") as handle:
      json.dump(data, handle, indent=2, allow_nan=False)
      handle.write("\n")
    os.replace(temporary, path)
    temporary = None
    state["startup"] = data["startup"]
    return True
  except Exception as error:
    info_dialog(f"Could not save trackpad startup settings.\n\n{error}")
    return False
  finally:
    if temporary is not None:
      try:
        os.unlink(temporary)
      except OSError:
        pass


def _trackpad_save_tick():
  _TRACKPAD_STATE["save_source"] = None
  if _TRACKPAD_STATE["startup"]:
    _save_trackpad_preferences()
  return False


def _schedule_trackpad_preferences_save():
  state = _TRACKPAD_STATE
  if state["startup"] and state["save_source"] is None:
    state["save_source"] = _coot_gui_repository_module("GLib").timeout_add(400, _trackpad_save_tick)


def _trackpad_startup_tick():
  state = _TRACKPAD_STATE
  if not state["startup"] or state["enabled"]:
    state["startup_source"] = None
    return False
  try:
    if _trackpad_model_areas():
      state["startup_source"] = None
      set_custom_trackpad_controls(True)
      return False
    # Wait only during startup, with a bound in case the drawing area never maps.
    state["startup_attempts"] += 1
    if state["startup_attempts"] < 20:
      return True
    print("WARNING:: trackpad startup skipped: no visible model drawing area")
  except Exception as error:
    print(f"WARNING:: could not start custom trackpad controls: {error}")
  state["startup_source"] = None
  return False


def _initialize_trackpad_preferences():
  state = _TRACKPAD_STATE
  try:
    with open(_trackpad_preferences_path(), encoding="utf-8") as handle:
      data = _trackpad_normalize_preferences(json.load(handle))
    for name in ("startup", "reverse", "sensitivity"):
      state[name] = data[name]
  except FileNotFoundError:
    return
  except Exception as error:
    print(f"WARNING:: could not read trackpad settings; using defaults: {error}")
    return
  if state["startup"] and state["startup_source"] is None:
    glib = _coot_gui_repository_module("GLib")
    if glib is not None:
      state["startup_attempts"] = 0
      state["startup_source"] = glib.timeout_add(250, _trackpad_startup_tick)


def _trackpad_model_areas():
  """Search only the main Coot view, not other windows such as Ramachandran plots."""
  gui = _coot_gui_module()
  if gui is None:
    return []
  pending = [gui.coot_gui_api.main_hbox()]
  areas = []
  while pending:
    widget = pending.pop()
    if isinstance(widget, Gtk.GLArea) and widget.get_mapped():
      areas.append(widget)
    child = widget.get_first_child()
    while child is not None:
      pending.append(child)
      child = child.get_next_sibling()
  if not areas:
    return []
  largest = max(area.get_width() * area.get_height() for area in areas)
  return [area for area in areas if area.get_width() * area.get_height() >= largest * 0.8]


def _trackpad_unit_vector(vector):
  length = math.sqrt(sum(float(v) ** 2 for v in vector))
  if not math.isfinite(length) or length < 1e-8:
    raise ValueError("Coot did not return usable screen vectors")
  return tuple(float(v) / length for v in vector)


def _trackpad_screen_frame():
  vectors = coot.screen_vectors_py()
  x = _trackpad_unit_vector(vectors[0])
  # Coot's screen Y vector points down, not up. Orthogonalize to limit drift.
  down = vectors[1]
  dot = sum(x[i] * float(down[i]) for i in range(3))
  y = _trackpad_unit_vector(tuple(-float(down[i]) + dot * x[i] for i in range(3)))
  z = (x[1]*y[2] - x[2]*y[1], x[2]*y[0] - x[0]*y[2], x[0]*y[1] - x[1]*y[0])
  return x, y, z


def _trackpad_frame_quaternion(frame):
  """Get XYZW from the view basis; older native quaternion getters return zero."""
  m = frame
  trace = m[0][0] + m[1][1] + m[2][2]
  if trace > 0:
    s = math.sqrt(trace + 1.0) * 2.0
    return ((m[2][1]-m[1][2])/s, (m[0][2]-m[2][0])/s, (m[1][0]-m[0][1])/s, s/4.0)
  i = max(range(3), key=lambda j: m[j][j])
  j, k = (i+1) % 3, (i+2) % 3
  s = math.sqrt(max(0.0, 1.0 + m[i][i] - m[j][j] - m[k][k])) * 2.0
  q = [0.0] * 4
  q[i], q[j], q[k], q[3] = s/4.0, (m[i][j]+m[j][i])/s, (m[i][k]+m[k][i])/s, (m[k][j]-m[j][k])/s
  return tuple(q)


def _trackpad_adjust_clipping(dy, sensitivity):
  if dy == 0:
    return
  front, back = float(coot.get_clipping_plane_front()), float(coot.get_clipping_plane_back())
  if not all(math.isfinite(v) for v in (front, back)):
    raise ValueError("Coot returned invalid clipping planes")
  perspective = bool(coot.use_perspective_projection_state())
  # Orthographic near/far planes use -0.1*zoom*front and +0.3*zoom*back.
  # Perspective getters instead return the actual near/far camera distances.
  a, b = (front, back) if perspective else (-0.1 * front, 0.3 * back)
  half = (b - a) / 2.0
  if half <= 0:
    return
  factor = math.exp(max(-2.0, min(2.0, dy * 0.01 * sensitivity)))
  if perspective:
    delta = half * (factor - 1.0)
    delta = max(-half * 0.99, min(delta, front - 2.01, 999.9 - back))
    target_front, target_back = front - delta, back + delta
  else:
    # Centre on the rotation centre, NOT the existing slab midpoint. Coot's
    # default front=back=1 has three times as much background as foreground.
    half = min(100.0, max(0.00003, half * factor))
    target_front, target_back = half / 0.1, half / 0.3
  coot.set_clipping_front(target_front)
  coot.set_clipping_back(target_back)
  if perspective:
    # Perspective also has eye-position limits not exposed by the API. If one
    # setter rejects a plane, restore both to preserve the requested geometry.
    actual_front = float(coot.get_clipping_plane_front())
    actual_back = float(coot.get_clipping_plane_back())
    if not math.isclose(front - actual_front, actual_back - back, rel_tol=1e-4, abs_tol=1e-4):
      coot.set_clipping_front(front)
      coot.set_clipping_back(back)


def _trackpad_adjust_proportional_radius(dy):
  if dy == 0:
    return
  state = _TRACKPAD_STATE
  step = _proportional_editing_radius_step_size()
  delta = state.get("radius_remainder", 0.0) - dy * 0.1
  steps = math.trunc(round(delta / step, 9))
  state["radius_remainder"] = delta - steps * step
  if steps:
    target = min(1000.0, max(0.0, PROPORTIONAL_EDITING_RADIUS + steps * step))
    _set_proportional_editing_radius_value(target)


def _trackpad_adjust_map_radius(dy):
  if dy == 0:
    return
  map_id = coot.scroll_wheel_map()
  if not coot.is_valid_map_molecule(map_id):
    return
  state = _TRACKPAD_STATE
  # Rebuilding map meshes is expensive; accumulate small swipes into 2 A steps.
  delta = state.get("map_radius_remainder", 0.0) - dy * 0.1
  steps = math.trunc(round(delta / 2.0, 9))
  state["map_radius_remainder"] = delta - steps * 2.0
  if steps:
    current = float(_current_map_radius(map_id))
    target = max(2.0, current + steps * 2.0)
    if target != current:
      _set_map_radius_both(target, map_id)
      add_status_bar_text(f"Map radius: {target:g} A")


def _trackpad_apply_motion(mode, value, area):
  edit = globals().get("_TRACKPAD_EDIT_STATE")
  if edit is not None and edit["move_selection"] and mode in ("rotate", "translate", "slab_translate"):
    _trackpad_edit_motion(edit, mode, value, area)
    return
  if mode == "radius":
    _trackpad_adjust_proportional_radius(value[1])
    return
  if mode == "map_radius":
    _trackpad_adjust_map_radius(value[1])
    return
  sensitivity = _TRACKPAD_STATE["sensitivity"]["clipping" if mode == "slab_translate" else mode]
  if mode == "zoom":
    coot.set_zoom(min(2000.0, max(0.2, float(value))))
    return
  dx, dy = value
  if mode == "clipping":
    _trackpad_adjust_clipping(dy, sensitivity)
    return
  if mode == "contour":
    map_id = coot.scroll_wheel_map()
    if coot.is_valid_map_molecule(map_id):
      level = coot.get_contour_level_in_sigma(map_id)
      coot.set_contour_level_in_sigma(map_id, max(0.01, level - dy * 0.01 * sensitivity))
    return
  frame = _trackpad_screen_frame()
  if mode in ("translate", "slab_translate"):
    scale = 0.6 * float(coot.zoom_factor()) / max(1, area.get_height()) * sensitivity
    centre = [coot.rotation_centre_position(i) for i in range(3)]
    qx, qy, qz, qw = _trackpad_frame_quaternion(frame)
    # set_rotation_centre animates moves and ignores displacements below 0.3 A.
    # The JSON setter jumps directly; its quaternion order is WXYZ, not XYZW.
    if mode == "slab_translate":
      # Moving the camera target in screen Z carries both world-space clip
      # planes and the rotation centre together. Plane settings stay unchanged.
      target = [centre[i] - scale * dy * frame[2][i] for i in range(3)]
    else:
      target = [centre[i] - scale * (dx * frame[0][i] - dy * frame[1][i]) for i in range(3)]
    coot.set_view_from_json(json.dumps({
      "rotation-centre": target,
      "quaternion": [qw, qx, qy, qz], "zoom": float(coot.zoom_factor())}))
    if coot.get_active_map_drag_flag():
      coot.update_maps()
    if coot.get_show_symmetry():
      coot.set_show_symmetry_master(1)
    return
  qx, qy, qz, qw = _trackpad_frame_quaternion(frame)
  ax, ay = math.radians(dy * 0.25 * sensitivity), math.radians(dx * 0.25 * sensitivity)
  angle = math.hypot(ax, ay)
  if angle < 1e-12:
    return
  factor = math.sin(angle / 2.0) / angle
  x, y, w = ax * factor, ay * factor, math.cos(angle / 2.0)
  q = (w*qx + x*qw + y*qz, w*qy + y*qw - x*qz,
       w*qz + x*qy - y*qx, w*qw - x*qx - y*qy)
  length = math.sqrt(sum(v*v for v in q))
  coot.set_view_quaternion(*(v / length for v in q))


def _trackpad_flush():
  state = _TRACKPAD_STATE
  state["source"] = None
  pending, state["pending"] = state["pending"], {}
  if not state["enabled"]:
    return False
  try:
    for (mode, area), value in pending.items():
      _trackpad_apply_motion(mode, value, area)
  except Exception as error:
    if globals().get("_TRACKPAD_EDIT_STATE") is not None:
      _trackpad_edit_finish(False)
    set_custom_trackpad_controls(False)
    add_status_bar_text(f"Trackpad controls disabled: {error}")
    traceback.print_exc()
  return False


def _trackpad_queue(mode, value, area):
  state = _TRACKPAD_STATE
  key = (mode, area)
  if mode == "zoom":
    state["pending"][key] = value
  else:
    previous = state["pending"].get(key, (0.0, 0.0))
    state["pending"][key] = tuple(a+b for a, b in zip(previous, value))
  if state["source"] is None:
    state["source"] = _coot_gui_repository_module("GLib").timeout_add(16, _trackpad_flush)


def _trackpad_cancel_pinch(area):
  _TRACKPAD_STATE["pinch"].pop(area, None)
  _TRACKPAD_STATE["pending"].pop(("zoom", area), None)


def _trackpad_scroll_begin(_controller, area):
  if _TRACKPAD_STATE["enabled"]:
    _TRACKPAD_STATE["scroll_active"][area] = True


def _trackpad_scroll_end(_controller, area):
  _TRACKPAD_STATE["scroll_active"].pop(area, None)


def _trackpad_claim_scroll(controller, area):
  _TRACKPAD_STATE["claimed_scroll"][area] = controller.get_current_event_time()
  return True


def _trackpad_native_scroll_guard(controller, dx, dy, area, record):
  """Never let a native wheel action repeat a captured trackpad event."""
  state = _TRACKPAD_STATE
  if not state["enabled"]:
    return False
  if state["claimed_scroll"].get(area) == controller.get_current_event_time():
    return True
  if _trackpad_scroll(controller, dx, dy, area):
    return True
  if not state["enabled"]:
    return False
  # Keep native wheel handling, including its real event/modifier context.
  # Only replay events the custom controls declined; block this guard to
  # avoid recursive emission and restore all handler states afterwards.
  controller.handler_block(record["gate"])
  try:
    for handler in record["handlers"]:
      controller.handler_unblock(handler)
    return bool(controller.emit("scroll", dx, dy))
  finally:
    for handler in record["handlers"]:
      controller.handler_block(handler)
    controller.handler_unblock(record["gate"])


def _trackpad_guard_native_scroll(controller, area):
  gobject = _coot_gui_repository_module("GObject")
  signal = gobject.signal_lookup("scroll", Gtk.EventControllerScroll)
  record = {"controller": controller, "handlers": [], "gate": None}
  _TRACKPAD_STATE["native_scroll"].append(record)
  mask = gobject.SignalMatchType.ID | gobject.SignalMatchType.UNBLOCKED
  while True:
    handler = gobject.signal_handler_find(controller, mask, signal, 0, None, None, None)
    if not handler:
      break
    if handler in record["handlers"]:
      raise RuntimeError("Coot's native wheel handler could not be suspended safely.")
    controller.handler_block(handler)
    record["handlers"].append(handler)
  record["gate"] = controller.connect("scroll", _trackpad_native_scroll_guard, area, record)


def _trackpad_scroll(controller, dx, dy, area):
  """Use scalar GTK signals: raw Gdk.Event arguments can arrive as None."""
  state = _TRACKPAD_STATE
  state["claimed_scroll"].pop(area, None)
  if not state["enabled"]:
    return False
  try:
    gdk = _coot_gui_repository_module("Gdk")
    modifiers = controller.get_current_event_state()
    if modifiers & (gdk.ModifierType.SHIFT_MASK | gdk.ModifierType.CONTROL_MASK |
                    gdk.ModifierType.ALT_MASK | gdk.ModifierType.META_MASK):
      _trackpad_cancel_pinch(area)
    buttons = gdk.ModifierType.BUTTON1_MASK | gdk.ModifierType.BUTTON2_MASK | gdk.ModifierType.BUTTON3_MASK
    if modifiers & (buttons | gdk.ModifierType.META_MASK):
      return False
    if (modifiers & gdk.ModifierType.SHIFT_MASK and modifiers & gdk.ModifierType.ALT_MASK
        and modifiers & gdk.ModifierType.CONTROL_MASK):
      return False
    dx, dy = float(dx), float(dy)
    # A zero-delta stop event must not reach Coot: native Shift-scroll treats
    # even dy=0 as a zoom step. Handle this before deciding wheel vs surface.
    if not math.isfinite(dx) or not math.isfinite(dy) or (dx == 0.0 and dy == 0.0):
      return _trackpad_claim_scroll(controller, area)
    # Device source is MOUSE even for a Mac trackpad; use scroll units instead.
    # Keep a continuous gesture captured even if GTK reports a wheel unit for
    # one sample; otherwise native Shift-scroll can insert a fixed zoom step.
    slab_translate = bool(modifiers & gdk.ModifierType.SHIFT_MASK and modifiers & gdk.ModifierType.CONTROL_MASK)
    radius = bool(modifiers & gdk.ModifierType.CONTROL_MASK and modifiers & gdk.ModifierType.ALT_MASK)
    map_radius = bool(modifiers & gdk.ModifierType.SHIFT_MASK and modifiers & gdk.ModifierType.ALT_MASK)
    # Reserve this chord even outside scroll-begin/end: wheel-unit samples
    # otherwise reach Coot's native Shift-Ctrl representation cycling.
    edit = globals().get("_TRACKPAD_EDIT_STATE")
    moving_selection = edit is not None and edit["move_selection"]
    # Some Mac swipe samples arrive as wheel events after a native dialog.
    # During selection movement, do not fall through to native contouring.
    if not (slab_translate or radius or map_radius or moving_selection) and controller.get_unit() != gdk.ScrollUnit.SURFACE and area not in state["scroll_active"]:
      return False
    direction = -1.0 if state["reverse"] else 1.0
    mode = ("map_radius" if map_radius else
            "radius" if radius else
            "slab_translate" if slab_translate else
            "contour" if modifiers & gdk.ModifierType.ALT_MASK else
            "translate" if modifiers & gdk.ModifierType.CONTROL_MASK else
            "clipping" if modifiers & gdk.ModifierType.SHIFT_MASK else "rotate")
    _trackpad_queue(mode, (dx * direction, dy * direction), area)
    return _trackpad_claim_scroll(controller, area)
  except Exception as error:
    set_custom_trackpad_controls(False)
    add_status_bar_text(f"Trackpad controls disabled: {error}")
    traceback.print_exc()
    return False


def _trackpad_pinch_begin(gesture, _sequence, area):
  if not _TRACKPAD_STATE["enabled"]:
    return
  try:
    gdk = _coot_gui_repository_module("Gdk")
    buttons = gdk.ModifierType.BUTTON1_MASK | gdk.ModifierType.BUTTON2_MASK | gdk.ModifierType.BUTTON3_MASK
    blocked = buttons | gdk.ModifierType.SHIFT_MASK | gdk.ModifierType.CONTROL_MASK | gdk.ModifierType.ALT_MASK | gdk.ModifierType.META_MASK
    if gesture.get_current_event_state() & blocked:
      _trackpad_cancel_pinch(area)
      return
    state = _TRACKPAD_STATE
    state["pinch"][area] = state["pending"].get(("zoom", area), float(coot.zoom_factor()))
  except Exception as error:
    set_custom_trackpad_controls(False)
    add_status_bar_text(f"Trackpad controls disabled: {error}")
    traceback.print_exc()


def _trackpad_pinch_scale(gesture, scale, area):
  state = _TRACKPAD_STATE
  if not state["enabled"] or area not in state["pinch"]:
    return
  try:
    gdk = _coot_gui_repository_module("Gdk")
    blocked = (gdk.ModifierType.SHIFT_MASK | gdk.ModifierType.CONTROL_MASK | gdk.ModifierType.ALT_MASK |
               gdk.ModifierType.META_MASK | gdk.ModifierType.BUTTON1_MASK |
               gdk.ModifierType.BUTTON2_MASK | gdk.ModifierType.BUTTON3_MASK)
    if gesture.get_current_event_state() & blocked:
      _trackpad_cancel_pinch(area)
      return
    scale = float(scale)
    if not math.isfinite(scale) or scale <= 0:
      return
    zoom = state["pinch"][area] * math.exp(max(-20.0, min(20.0, -state["sensitivity"]["zoom"] * math.log(scale))))
    _trackpad_queue("zoom", zoom, area)
  except Exception as error:
    set_custom_trackpad_controls(False)
    add_status_bar_text(f"Trackpad controls disabled: {error}")
    traceback.print_exc()


def _trackpad_pinch_end(_gesture, _sequence, area):
  _TRACKPAD_STATE["pinch"].pop(area, None)


def set_custom_trackpad_controls(enabled=True):
  """Install capture handlers and temporarily suspend Coot's pinch recognizers."""
  state = _TRACKPAD_STATE
  if not enabled:
    state["enabled"] = False
    if state["startup_source"] is not None:
      _coot_gui_repository_module("GLib").source_remove(state["startup_source"])
      state["startup_source"] = None
    if state["source"] is not None:
      _coot_gui_repository_module("GLib").source_remove(state["source"])
      state["source"] = None
    state["pending"].clear()
    state["radius_remainder"] = 0.0
    state["map_radius_remainder"] = 0.0
    state["pinch"].clear()
    state["scroll_active"].clear()
    state["claimed_scroll"].clear()
    while state["bindings"]:
      area, controller = state["bindings"].pop()
      try:
        area.remove_controller(controller)
      except Exception as error:
        print("Could not detach trackpad controller:", error)
    while state["native_scroll"]:
      record = state["native_scroll"].pop()
      controller = record["controller"]
      if record["gate"] is not None:
        controller.disconnect(record["gate"])
      for handler in record["handlers"]:
        controller.handler_unblock(handler)
    while state["native_zoom"]:
      controller, phase = state["native_zoom"].pop()
      controller.set_propagation_phase(phase)
    window = state.get("window")
    if window is not None and hasattr(window, "trackpad_widgets"):
      checkbox = window.trackpad_widgets["enabled"]
      if checkbox.get_active():
        checkbox.set_active(False)
    add_status_bar_text("Native Coot trackpad controls restored")
    return True
  if state["enabled"]:
    return True
  try:
    gdk = _coot_gui_repository_module("Gdk")
    if (gdk is None or not hasattr(gdk, "ScrollUnit") or
        not hasattr(Gtk.EventControllerScroll, "get_unit")):
      raise RuntimeError("This mode needs GTK 4.8 or newer with Python event support")
    if not callable(getattr(coot, "set_view_from_json", None)):
      raise RuntimeError("This Coot build does not support direct camera panning")
    areas = _trackpad_model_areas()
    if not areas:
      raise RuntimeError("No visible Coot model drawing area was found")
    for area in areas:
      controllers = area.observe_controllers()
      for index in range(controllers.get_n_items()):
        native = controllers.get_item(index)
        if isinstance(native, Gtk.GestureZoom):
          state["native_zoom"].append((native, native.get_propagation_phase()))
          native.set_propagation_phase(Gtk.PropagationPhase.NONE)
        elif isinstance(native, Gtk.EventControllerScroll):
          _trackpad_guard_native_scroll(native, area)
      controller = Gtk.EventControllerScroll.new(Gtk.EventControllerScrollFlags.BOTH_AXES)
      controller.set_propagation_phase(Gtk.PropagationPhase.CAPTURE)
      controller.connect("scroll", _trackpad_scroll, area)
      controller.connect("scroll-begin", _trackpad_scroll_begin, area)
      controller.connect("scroll-end", _trackpad_scroll_end, area)
      area.add_controller(controller)
      state["bindings"].append((area, controller))
      pinch = Gtk.GestureZoom.new()
      pinch.set_propagation_phase(Gtk.PropagationPhase.CAPTURE)
      pinch.connect("begin", _trackpad_pinch_begin, area)
      pinch.connect("scale-changed", _trackpad_pinch_scale, area)
      pinch.connect("end", _trackpad_pinch_end, area)
      pinch.connect("cancel", _trackpad_pinch_end, area)
      area.add_controller(pinch)
      state["bindings"].append((area, pinch))
    state["enabled"] = True
    add_status_bar_text("Trackpad: pinch zoom; swipe rotate; Ctrl pan; Option contour; Ctrl-Option editing radius; Shift slab width; Shift-Ctrl slab Z")
    return True
  except Exception as error:
    set_custom_trackpad_controls(False)
    info_dialog(f"Could not enable custom trackpad controls.\n\n{error}")
    return False


def _trackpad_edit_clear_pending():
  state = _TRACKPAD_STATE
  if state["source"] is not None:
    _coot_gui_repository_module("GLib").source_remove(state["source"])
    state["source"] = None
  state["pending"].clear()
  state["pinch"].clear()


def _trackpad_native_widgets():
  """Find native controls once per preview, never in the gesture hot path."""
  roots = Gtk.Window.get_toplevels()
  pending = [roots.get_item(i) for i in range(roots.get_n_items())]
  found = {}
  while pending:
    widget = pending.pop()
    name = Gtk.Buildable.get_buildable_id(widget)
    if name:
      found[name] = widget
    child = widget.get_first_child()
    while child is not None:
      pending.append(child)
      child = child.get_next_sibling()
  return found


def _trackpad_native_close_requested(dialog):
  """Keep Coot's reusable dialog alive and cancel its abandoned preview."""
  edit = globals().get("_TRACKPAD_EDIT_STATE")
  if edit is not None and edit.get("owns_native") and edit.get("native") is dialog:
    _trackpad_edit_reset()
  else:
    coot.clear_up_moving_atoms()
    coot.rot_trans_reset_previous()
  dialog.set_visible(False)
  coot.normal_cursor()
  return True


def _install_trackpad_native_close_guard(widgets=None):
  widgets = _trackpad_native_widgets() if widgets is None else widgets
  dialog = widgets.get("rotate_translate_obj_dialog")
  if dialog is not None and not getattr(dialog, "coot_trimmings_close_guard", None):
    dialog.coot_trimmings_close_guard = dialog.connect("close-request", _trackpad_native_close_requested)


def _trackpad_edit_update_buttons(edit):
  button = edit.get("rigid_fit_button")
  if button is not None:
    button.set_sensitive(bool(edit.get("idle")) and not edit.get("fit_busy", False))


def _trackpad_edit_begin_native(edit):
  """Start Coot's own moving-atoms buffer, without a new model molecule."""
  if coot.refinement_already_ongoing_p():
    raise RuntimeError("Accept or cancel the current intermediate atoms first.")
  action = _coot_gui_module().coot_gui_api.application().lookup_action("rotate_translate_" + edit["scope"])
  if action is None:
    raise RuntimeError("This Coot build does not expose the native Rotate/Translate action.")
  # Native actions target the nearest atom. Temporarily put the stored target
  # at the view centre and isolate picking, then restore everything in this UI turn.
  centre = [coot.rotation_centre_position(i) for i in range(3)]
  qx, qy, qz, qw = _trackpad_frame_quaternion(_trackpad_screen_frame())
  view = {"rotation-centre": centre, "quaternion": [qw, qx, qy, qz], "zoom": float(coot.zoom_factor())}
  picking = [(m, coot.mol_is_active(m), coot.mol_is_displayed(m)) for m in model_molecule_list()]
  try:
    if edit["scope"] != "residue_range":
      atoms = coot.residue_info_py(edit["source"], *edit["specs"][0]) or []
      if not atoms:
        raise RuntimeError("The selected residue is no longer available.")
      atom = next((a for a in atoms if not a[0][1] and a[0][0].strip() == "CA"), atoms[0])
      for model, _active, _displayed in picking:
        coot.set_mol_active(model, int(model == edit["source"]))
      coot.set_mol_displayed(edit["source"], 1)
      coot.set_view_from_json(json.dumps(dict(view, **{"rotation-centre": atom[2]})))
    elif not edit.get("range_valid", False):
      raise RuntimeError("The native residue range changed. Close this panel and select it again.")
    coot.rot_trans_reset_previous()
    # Coot exposes a setter but no getter for this native rotation preference.
    coot.set_rotate_translate_zone_rotates_about_zone_centre(1)
    action.activate(None)
    edit["owns_native"] = True
  finally:
    coot.set_view_from_json(json.dumps(view))
    for model, active, displayed in picking:
      if valid_model_molecule_qm(model):
        coot.set_mol_active(model, active)
        coot.set_mol_displayed(model, displayed)
  widgets = _trackpad_native_widgets()
  _install_trackpad_native_close_guard(widgets)
  dialog = widgets.get("rotate_translate_obj_dialog")
  edit["native"] = dialog
  if dialog is None or not dialog.get_visible() or not coot.refinement_already_ongoing_p():
    _trackpad_edit_discard_native(edit)
    raise RuntimeError("Coot could not create the native movement preview.")
  edit["adjustments"] = [widgets["rotate_translate_obj_" + name + "_hscale"].get_adjustment()
                         for name in ("xtrans", "ytrans", "ztrans", "xrot", "yrot", "zrot")]
  dialog.set_visible(False)
  edit["original"] = _atom_snapshot_for_specs(edit["source"], edit["specs"])
  edit["dirty"] = False
  edit["idle"] = False
  _trackpad_edit_update_buttons(edit)


def _trackpad_edit_discard_native(edit):
  if edit.get("owns_native"):
    coot.clear_up_moving_atoms()
    coot.rot_trans_reset_previous()
    if edit.get("native") is not None:
      edit["native"].set_visible(False)
  edit["owns_native"] = False
  edit["adjustments"] = []
  edit["idle"] = True
  _trackpad_edit_update_buttons(edit)


def _trackpad_edit_check(edit, check_coordinates=True):
  if not valid_model_molecule_qm(edit["source"]):
    raise RuntimeError("The selected model was closed.")
  if not edit.get("adjustments") or not coot.refinement_already_ongoing_p():
    raise RuntimeError("The native movement preview ended. Close this panel and start again.")
  if check_coordinates and _atom_snapshot_for_specs(edit["source"], edit["specs"]) != edit["original"]:
    raise RuntimeError("The selected atoms changed outside this panel. Close and reopen it.")


def _trackpad_edit_adjust(edit, increments):
  adjustments = edit["adjustments"]
  # Native sliders are incremental and bounded. Rebase their remembered
  # values before reaching a limit; sentinel + zero means no physical move.
  steps = max(1, int(math.ceil(max(abs(v) for v in increments) / 90.0)))
  for _ in range(steps):
    if any(not a.get_lower() < a.get_value() + v/steps < a.get_upper()
           for a, v in zip(adjustments, increments)):
      coot.rot_trans_reset_previous()
      for adjustment in adjustments:
        adjustment.set_value(0.0)
    for adjustment, value in zip(adjustments, increments):
      if value:
        adjustment.set_value(adjustment.get_value() + value/steps)
  edit["dirty"] = True


def _trackpad_edit_motion(edit, mode, value, area):
  if edit.get("idle"):
    if not valid_model_molecule_qm(edit["source"]):
      raise RuntimeError("The selected model was closed.")
    _trackpad_edit_begin_native(edit)
  _trackpad_edit_check(edit, check_coordinates=False)
  dx, dy = value
  if not all(math.isfinite(v) for v in (dx, dy)) or max(abs(dx), abs(dy)) > 100000:
    raise ValueError("Invalid trackpad movement.")
  key = "clipping" if mode == "slab_translate" else mode
  sensitivity = _TRACKPAD_STATE["sensitivity"][key]
  increments = [0.0]*6
  if mode == "rotate":
    scale = math.radians(0.25) * sensitivity / 0.018
    increments[3], increments[4] = dy * scale, -dx * scale
  else:
    scale = 300.0 / max(1, area.get_height()) * sensitivity
    if mode == "slab_translate":
      increments[2] = dy * scale
    else:
      increments[0], increments[1] = dx * scale, -dy * scale
  _trackpad_edit_adjust(edit, increments)


def _trackpad_edit_reset():
  edit = _TRACKPAD_EDIT_STATE
  if edit is None:
    return
  _trackpad_edit_clear_pending()
  if edit.get("idle"):
    add_status_bar_text("No unapplied movement to reset.")
    return
  try:
    _trackpad_edit_check(edit)
    _trackpad_edit_discard_native(edit)
    add_status_bar_text("Movement reset to the last applied coordinates.")
  except Exception as error:
    _trackpad_edit_finish(False)
    info_dialog(str(error))


def _trackpad_edit_fit_selection(edit):
  """Rigid fitting accepts one MMDB CID, not Coot's '||' selection extension."""
  scope, specs = edit["scope"], edit["specs"]
  if scope == "molecule":
    return "/"
  if not specs or len({spec[0] for spec in specs}) != 1:
    raise ValueError("The stored fitting selection must be in one chain.")
  chain = specs[0][0]
  if not chain or any(c in chain for c in "/|*,[]"):
    raise ValueError("This chain identifier cannot be represented safely by the rigid-fit API.")
  if scope == "chain":
    return "//" + chain
  if scope == "residue" and len(specs) == 1:
    return _residue_spec_to_cid_compat(specs[0])
  if scope == "residue_range":
    start, end = min(s[1] for s in specs), max(s[1] for s in specs)
    actual = {tuple(s) for s in _tq_chain_residue_specs(edit["source"], chain) if start <= s[1] <= end}
    if actual != set(map(tuple, specs)):
      raise ValueError("The stored range no longer matches a complete native residue range. Select it again.")
    return f"//{chain}/{start}-{end}"
  raise ValueError("Unrecognised rigid-fit selection scope.")


def _trackpad_edit_rigid_fit():
  edit = _TRACKPAD_EDIT_STATE
  if edit is None or edit.get("fit_busy"):
    return False
  pending_movement = edit["move_selection"] and any(
    mode in ("rotate", "translate", "slab_translate") for mode, _area in _TRACKPAD_STATE["pending"])
  if not edit.get("idle") or pending_movement:
    add_status_bar_text("Apply or reset the manual movement before rigid fitting.")
    return False
  if coot.refinement_already_ongoing_p():
    add_status_bar_text("Accept or cancel the current intermediate atoms before rigid fitting.")
    return False
  previous_immediate = None
  started_fit = False
  try:
    if not valid_model_molecule_qm(edit["source"]):
      raise RuntimeError("The selected model was closed.")
    fit = getattr(coot, "rigid_body_refine_by_atom_selection", None)
    if not callable(fit):
      raise RuntimeError("This Coot build cannot rigid-fit an explicit atom selection.")
    map_id = coot.imol_refinement_map()
    if not coot.is_valid_map_molecule(map_id):
      raise RuntimeError("Select a refinement map before rigid fitting.")
    if coot.backup_state(edit["source"]) != 1:
      raise RuntimeError("Enable Coot backups for this model so the fit can be undone.")
    _atom_snapshot_for_specs(edit["source"], edit["specs"])
    selection = _trackpad_edit_fit_selection(edit)
    _trackpad_edit_clear_pending()
    edit["fit_busy"] = True
    _trackpad_edit_update_buttons(edit)
    previous_immediate = coot.refinement_immediate_replacement_state()
    coot.set_refinement_immediate_replacement(1)
    started_fit = True
    # This API creates fitted intermediate atoms; it does not commit them.
    fit(edit["source"], selection)
    if not coot.refinement_already_ongoing_p():
      raise RuntimeError("Coot did not produce a rigid-fit solution; no fit was applied.")
    with _grouped_model_edit(edit["source"]):
      coot.accept_regularizement()
    add_status_bar_text(f"Rigid fit applied to the stored selection using map #{map_id}. Swipe to continue; Undo restores the fit.")
    return True
  except Exception as error:
    info_dialog(f"Could not rigid-fit the selection.\n\n{error}")
    return False
  finally:
    try:
      if started_fit:
        coot.clear_up_moving_atoms()
    finally:
      if previous_immediate is not None:
        coot.set_refinement_immediate_replacement(previous_immediate)
      edit["fit_busy"] = False
      _trackpad_edit_update_buttons(edit)


def _trackpad_edit_finish(accept=False, keep_open=False):
  global _TRACKPAD_EDIT_STATE
  edit = _TRACKPAD_EDIT_STATE
  if edit is None:
    return False
  if accept and _TRACKPAD_STATE["pending"]:
    source = _TRACKPAD_STATE["source"]
    if source is not None:
      _coot_gui_repository_module("GLib").source_remove(source)
    _trackpad_flush()
    if _TRACKPAD_EDIT_STATE is not edit:
      return False
  else:
    _trackpad_edit_clear_pending()
  if accept:
    try:
      if not edit.get("idle"):
        _trackpad_edit_check(edit)
        if edit["dirty"]:
          with _grouped_model_edit(edit["source"]):
            # Calling the Python API retains the GIL. Emitting the native OK
            # button through GI can release it while Coot runs Python hooks.
            coot.accept_regularizement()
        _trackpad_edit_discard_native(edit)
      if keep_open:
        add_status_bar_text("Movement applied. Rigid fit is available; swipe to resume moving the same selection.")
        return True
    except Exception as error:
      info_dialog(f"Could not complete trackpad movement.\n\n{error}\nAny applied edits remain available in Coot Undo.")
      return False
  _TRACKPAD_EDIT_STATE = None
  try:
    _trackpad_edit_discard_native(edit)
  finally:
    if edit.get("range_watch") is not None:
      widget, handler = edit["range_watch"]
      widget.disconnect(handler)
    if edit.get("window") is not None:
      edit["window"].destroy()
    if not edit["was_enabled"]:
      set_custom_trackpad_controls(False)
  add_status_bar_text("Trackpad movement applied; one Undo restores it." if accept else
                      "Trackpad panel closed; unapplied movement discarded.")
  return True


def _trackpad_edit_panel(edit):
  window = Gtk.Window(title="Trackpad Rotate/Translate")
  edit["window"] = window
  parent = _coot_gui_module().coot_gui_api.main_toolbar().get_root()
  if isinstance(parent, Gtk.Window):
    window.set_transient_for(parent)
  box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=10)
  for edge in ("start", "end", "top", "bottom"):
    getattr(box, "set_margin_" + edge)(12)
  label = Gtk.Label(label=edit["label"] + "\nPivot: selection centre")
  label.set_xalign(0)
  box.append(label)
  instructions = Gtk.Label(label="Swipe: rotate | Ctrl: X/Y | Shift-Ctrl: Z\nPinch: zoom | Option: contour | Ctrl-Option: editing radius\nApply, then rigid fit or swipe again. Close discards unapplied moves.")
  instructions.set_xalign(0)
  box.append(instructions)
  view = Gtk.CheckButton(label="Move view instead of selection")
  def change_view(widget):
    _trackpad_edit_clear_pending()
    edit["move_selection"] = not widget.get_active()
  view.connect("toggled", change_view)
  box.append(view)
  buttons = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
  for name, callback in (("Rigid-fit selection", lambda *_: _trackpad_edit_rigid_fit()),
                         ("Reset", lambda *_: _trackpad_edit_reset()),
                         ("Apply", lambda *_: _trackpad_edit_finish(True, keep_open=True)),
                         ("Close", lambda *_: _trackpad_edit_finish(False))):
    button = Gtk.Button(label=name)
    button.connect("clicked", callback)
    if name == "Rigid-fit selection":
      edit["rigid_fit_button"] = button
      button.set_tooltip_text("Fit the stored selection to the selected refinement map. Apply or reset manual movement first.")
    buttons.append(button)
  box.append(buttons)
  window.set_child(box)
  window.connect("close-request", lambda *_: _trackpad_edit_finish(False) or True)
  window.trackpad_edit_widgets = {"view": view, "buttons": buttons, "rigid_fit": edit["rigid_fit_button"]}
  _trackpad_edit_update_buttons(edit)
  window.present()
  return window


def _trackpad_pick_native_range():
  """Let Coot record the range itself, observing real picks for our markers."""
  if coot.refinement_already_ongoing_p():
    add_status_bar_text("Accept or cancel the current intermediate atoms first.")
    return None
  _cancel_script_clicks()
  _clear_annotation_range_pick()
  coot.clear_pending_picks()
  widgets = _trackpad_native_widgets()
  toggle = widgets.get("select_residue_range_toggle_button")
  if toggle is None or not callable(getattr(coot, "select_atom_under_pointer_py", None)):
    info_dialog("This Coot build cannot expose native range picking for trackpad movement.")
    return None
  token, picks, markers, handlers = object(), [], [], []
  def cleanup():
    for widget, handler in handlers:
      widget.disconnect(handler)
    handlers.clear()
    toggle.set_active(False)
    while markers:
      coot.close_generic_object(markers.pop())
    coot.graphics_draw()
  def clicked(controller, n_press, _x, _y):
    if _SCRIPT_CLICK_STATE.get("token") is not token or n_press != 1 or controller.get_current_button() != 1:
      return
    # This handler runs AFTER Coot has picked/stored the endpoint. Unlike the
    # user-defined-click callback, this API returns an explicit molecule ID.
    selected = coot.select_atom_under_pointer_py()
    if not selected:
      return
    model, atom = selected
    spec = [model, *atom[-5:]]
    try:
      spec = _resolve_script_click(spec)
      picks.append(spec)
      _mark_register_endpoint(spec, "cyan" if len(picks) == 1 else "orange", markers, "Trackpad range endpoint")
      add_status_bar_text(_navigation_status_bar_label(model, _click_spec_chain_id(spec),
        _click_spec_res_no(spec), _click_spec_ins_code(spec)).replace("Current residue:", "Clicked residue:", 1)
        + (" | Click the end residue." if len(picks) == 1 else " | Range selected."))
      if len(picks) == 2:
        specs = _clicked_residue_specs(*picks)
        if not specs:
          _cancel_script_clicks()
          return
        native_specs = [tuple(s) for s in _tq_chain_residue_specs(model, specs[0][0])
                        if specs[0][1] <= s[1] <= specs[-1][1]]
        # Native range selection includes all insertion codes at its endpoints.
        if set(native_specs) != set(map(tuple, specs)):
          raise ValueError("Native range movement cannot reliably isolate insertion-coded endpoints; use residue or chain movement.")
        _SCRIPT_CLICK_STATE.clear()
        cleanup()
        if specs:
          _start_trackpad_edit(model, specs, f"#{model} | {specs[0][0]}:{specs[0][1]}-{specs[-1][1]}",
                               scope="residue_range", range_valid=True)
    except Exception as error:
      _cancel_script_clicks()
      info_dialog(f"Could not select a movement range.\n\n{error}")
  try:
    for area in _trackpad_model_areas():
      controllers = area.observe_controllers()
      for i in range(controllers.get_n_items()):
        controller = controllers.get_item(i)
        if isinstance(controller, Gtk.GestureClick) and controller.get_button() in (0, 1):
          handlers.append((controller, controller.connect_after("pressed", clicked)))
    if not handlers:
      raise RuntimeError("Could not find Coot's native atom-picking controller.")
    toggle.set_active(False)
    toggle.set_active(True)
    _SCRIPT_CLICK_STATE.update(token=token, picking=False, cleanup=cleanup)
    add_status_bar_text("Trackpad movement: click the start and end residues. Tab cancels.")
    return token
  except Exception as error:
    cleanup()
    info_dialog(str(error))
    return None


def _start_trackpad_edit(mol_id, specs, label, scope="residue", range_valid=False):
  global _TRACKPAD_EDIT_STATE
  if _TRACKPAD_EDIT_STATE is not None:
    _TRACKPAD_EDIT_STATE["window"].present()
    add_status_bar_text("Close the current movement panel before choosing a different selection.")
    return None
  edit = None
  was_enabled = _TRACKPAD_STATE["enabled"]
  try:
    if coot.refinement_already_ongoing_p():
      raise RuntimeError("Accept or cancel the current refinement/Rotate-Translate operation first.")
    if coot.backup_state(mol_id) != 1:
      raise RuntimeError("Enable Coot backups for this model so Apply can be undone.")
    specs = [tuple(spec) for spec in specs]
    original = _atom_snapshot_for_specs(mol_id, specs)
    if not set_custom_trackpad_controls(True):
      raise RuntimeError("Custom trackpad controls could not be enabled.")
    _trackpad_edit_clear_pending()
    edit = {"source": mol_id, "specs": specs, "scope": scope, "range_valid": range_valid,
            "move_selection": True, "adjustments": [], "was_enabled": was_enabled, "window": None, "label": label,
            "idle": True, "dirty": False, "owns_native": False, "original": original, "fit_busy": False}
    _TRACKPAD_EDIT_STATE = edit
    if scope == "residue_range":
      toggle = _trackpad_native_widgets()["select_residue_range_toggle_button"]
      def range_changed(widget):
        if widget.get_active():
          edit["range_valid"] = False
      edit["range_watch"] = (toggle, toggle.connect("toggled", range_changed))
    _trackpad_edit_panel(edit)
    add_status_bar_text("Move selection: swipe rotates; Ctrl moves X/Y; Shift-Ctrl moves Z. Apply keeps the panel open.")
    return edit
  except Exception as error:
    if edit is not None:
      _trackpad_edit_finish(False)
    else:
      if not was_enabled and _TRACKPAD_STATE["enabled"]:
        set_custom_trackpad_controls(False)
    info_dialog(f"Trackpad movement could not start.\n\n{error}")
    return None


def trackpad_rotate_translate(scope="residue"):
  if _TRACKPAD_EDIT_STATE is not None:
    _TRACKPAD_EDIT_STATE["window"].present()
    add_status_bar_text("Close the current movement panel before choosing a different selection.")
    return None
  if scope == "range":
    return _trackpad_pick_native_range()
  active = _active_residue_or_status()
  if not active:
    return None
  mol_id, chain, number, ins = active[:4]
  if scope == "residue":
    specs, label = [(chain, number, ins)], f"#{mol_id} | Residue {chain}:{number}{ins}"
  elif scope == "chain":
    specs, label = _tq_chain_residue_specs(mol_id, chain), f"#{mol_id} | Chain {chain}"
  elif scope == "molecule":
    specs, label = _all_residue_specs_for_colouring(mol_id), f"Molecule #{mol_id}"
  else:
    raise ValueError("Choose residue, range, chain or molecule.")
  return _start_trackpad_edit(mol_id, specs, label, scope=scope)


def custom_trackpad_controls_dialog():
  """Adjust current navigation and optionally remember it for startup."""
  state = _TRACKPAD_STATE
  if state["window"] is not None:
    state["window"].present()
    return state["window"]
  window = Gtk.Window(title="Trackpad controls")
  state["window"] = window
  window.set_default_size(460, -1)
  gui = _coot_gui_module()
  parent = gui.coot_gui_api.main_toolbar().get_root() if gui else None
  if isinstance(parent, Gtk.Window):
    window.set_transient_for(parent)
  box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=10)
  for edge in ("start", "end", "top", "bottom"):
    getattr(box, "set_margin_" + edge)(12)
  enabled = Gtk.CheckButton(label="Enable custom trackpad controls (this session)")
  enabled.set_active(state["enabled"])
  box.append(enabled)
  startup = Gtk.CheckButton(label="Enable automatically on startup")
  startup.set_active(state["startup"])
  box.append(startup)
  label = Gtk.Label(label="Pinch: zoom   |   Swipe: rotate\nCtrl + swipe: translate   |   Option + swipe: contour\nCtrl + Option + swipe: proportional editing radius\nShift + Option + swipe: map radius (2 A steps)\nShift + swipe: slab width   |   Shift + Ctrl + swipe: slab Z")
  label.set_xalign(0)
  box.append(label)
  grid = Gtk.Grid(column_spacing=12, row_spacing=6)
  sliders = {}
  for row, (mode, name) in enumerate((("zoom", "Zoom"), ("rotate", "Rotation"),
                                     ("translate", "Translation"), ("contour", "Contour"),
                                     ("clipping", "Clipping (width/Z)"))):
    grid.attach(Gtk.Label(label=name), 0, row, 1, 1)
    slider = Gtk.Scale.new_with_range(Gtk.Orientation.HORIZONTAL, 0.1, 3.0, 0.1)
    slider.set_digits(1)
    slider.set_draw_value(True)
    slider.set_value_pos(Gtk.PositionType.RIGHT)
    slider.set_hexpand(True)
    slider.set_value(state["sensitivity"][mode])
    def sensitivity_changed(widget, key=mode):
      state["sensitivity"][key] = widget.get_value()
      _schedule_trackpad_preferences_save()
    slider.connect("value-changed", sensitivity_changed)
    grid.attach(slider, 1, row, 1, 1)
    sliders[mode] = slider
  box.append(grid)
  reverse = Gtk.CheckButton(label="Reverse swipe direction")
  reverse.set_active(state["reverse"])
  def reverse_changed(widget):
    state["reverse"] = widget.get_active()
    _schedule_trackpad_preferences_save()
  reverse.connect("toggled", reverse_changed)
  box.append(reverse)
  label = Gtk.Label(label="Sensitivity multipliers; 2.0 = default.\nMouse clicks and ordinary wheel input are unchanged.")
  label.set_xalign(0)
  box.append(label)
  buttons = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
  restore = Gtk.Button(label="Restore native controls")
  close = Gtk.Button(label="Close")
  buttons.append(restore)
  buttons.append(close)
  box.append(buttons)
  def toggle(widget):
    if not set_custom_trackpad_controls(widget.get_active()):
      widget.set_active(False)
  enabled.connect("toggled", toggle)
  def startup_changed(widget):
    wanted = widget.get_active()
    if wanted != state["startup"] and not _save_trackpad_preferences(wanted):
      widget.set_active(state["startup"])
  startup.connect("toggled", startup_changed)
  restore.connect("clicked", lambda *_: enabled.set_active(False))
  close.connect("clicked", lambda *_: window.close())
  def closed(*_args):
    if state["save_source"] is not None:
      _coot_gui_repository_module("GLib").source_remove(state["save_source"])
      _trackpad_save_tick()
    state["window"] = None
    return False
  window.connect("close-request", closed)
  window.trackpad_widgets = {"enabled": enabled, "startup": startup, "sliders": sliders, "reverse": reverse,
                            "restore": restore, "close": close}
  window.set_child(box)
  window.present()
  return window


def _set_native_high_contrast_model_defaults(show_dialog=True):
  """Apply Coot's native defaults for subsequently created models."""
  names = ("set_default_model_material_ambient", "set_default_model_material_diffuse",
           "set_default_model_material_specular", "set_default_bond_thickness",
           "set_use_variable_bond_thickness", "set_do_GL_lighting")
  missing = [name for name in names if not callable(getattr(coot, name, None))]
  if missing:
    if show_dialog:
      info_dialog("This Coot build cannot set the high-contrast model defaults.\\n\\n"
                  "Missing API: " + ", ".join(missing))
    return False
  try:
    coot.set_default_model_material_ambient(*MODEL_AMBIENT_ONLY_AMBIENT)
    coot.set_default_model_material_diffuse(*MODEL_AMBIENT_ONLY_DIFFUSE)
    # Coot 1.3.4's default-material API takes the molecule-type flag first.
    coot.set_default_model_material_specular(0, *MODEL_AMBIENT_ONLY_SPECULAR)
    coot.set_default_bond_thickness(HIGH_CONTRAST_BOND_THICKNESS)
    coot.set_use_variable_bond_thickness(0)
    coot.set_do_GL_lighting(1)
  except Exception as error:
    if show_dialog:
      info_dialog(f"Could not finish setting model defaults; some settings may have changed.\\n\\n{error}")
    else:
      print(f"WARNING:: could not restore high-contrast model defaults: {error}")
    return False
  return True


def default_models_to_high_contrast_display():
  """Persist and apply high-contrast defaults for subsequently created models."""
  if not _set_native_high_contrast_model_defaults(show_dialog=True):
    return False
  if not _save_display_default("high_contrast_models", True):
    return False
  message = "New models will default to high-contrast display, including future Coot sessions."
  add_status_bar_text(message)
  info_dialog(message + "\\n\\nExisting model materials are unchanged. "
              "Global lighting is enabled and variable bond thickness is disabled.")
  return True


def _set_native_ambient_map_defaults(show_dialog=True):
  """Apply the ambient-map defaults natively exposed by Coot 1.3.4."""
  # The uploaded Coot 1.3.4 coot.py exposes exactly these two map-default
  # material setters.  Specular, Fresnel, and colour remain per-map settings.
  required = (
    "set_default_map_material_ambient",
    "set_default_map_material_diffuse",
  )
  missing = [name for name in required if not callable(getattr(coot, name, None))]
  if missing:
    if show_dialog:
      info_dialog("This Coot build cannot set the ambient map material defaults.\n\n"
                  "Missing API: " + ", ".join(missing) +
                  "\nThese native default-material functions require Coot 1.3.4 or newer.")
    return False
  try:
    coot.set_default_map_material_ambient(*EM_MAP_AMBIENT)
    coot.set_default_map_material_diffuse(*EM_MAP_DIFFUSE)
  except Exception as error:
    if show_dialog:
      info_dialog(f"Could not finish setting map defaults; some settings may have changed.\n\n{error}")
    else:
      print(f"WARNING:: could not restore ambient map defaults: {error}")
    return False
  return True


def default_maps_to_ambient_display():
  """Persist and apply ambient defaults for subsequently created maps."""
  if not _set_native_ambient_map_defaults(show_dialog=True):
    return False
  if not _save_display_default("ambient_maps", True):
    return False
  message = "New maps will default to ambient material display, including future Coot sessions."
  add_status_bar_text(message)
  info_dialog(
    message +
    "\n\nCoot 1.3.4 exposes native map defaults for ambient and diffuse material only. "
    "Specular = 0, Fresnel off, and dark blue (#1a5fb4) are therefore still applied "
    "by the Restyle/resample map command rather than as native defaults for every new map."
    "\n\nExisting maps are unchanged."
  )
  return True


def _apply_persistent_display_defaults():
  """Restore startup values that Coot may overwrite during initialization.

  Coot can restore parts of its own display state after startup scripts have
  been evaluated, so this function is also scheduled once on the GTK event
  loop below.  Map radius is reapplied unconditionally because it is a normal
  startup preference, not part of either optional display-style preset.
  Returning False makes it suitable as a GLib timeout callback.
  """
  # Coot's later state/startup restoration can reset these to its built-in
  # radius (commonly 10 or 20 A depending on context).  Reassert both here.
  set_map_radius(STARTUP_MAP_RADIUS)
  set_map_radius_em(STARTUP_MAP_RADIUS_EM)

  defaults = _load_display_defaults()
  if defaults.get("high_contrast_models"):
    _set_native_high_contrast_model_defaults(show_dialog=False)
    # Models supplied on Coot's command line can be instantiated before this
    # startup script runs, so they may already have copied Coot's built-in
    # bond thickness (typically 5) even though the global default is now 2.
    # This callback is one-shot: correct models that already exist at startup
    # without installing a recurring molecule watcher.
    try:
      for imol in model_molecule_list():
        set_bond_thickness(imol, HIGH_CONTRAST_BOND_THICKNESS)
    except Exception as error:
      print(f"WARNING:: could not apply startup bond thickness to existing models: {error}")
  if defaults.get("ambient_maps"):
    _set_native_ambient_map_defaults(show_dialog=False)
  if defaults.get("black_background"):
    _set_black_background(show_dialog=False)
  _set_refinement_ramachandran_overlay(defaults.get("refinement_ramachandran_overlay", True), show_dialog=False)
  # Do not reapply contact-dot preferences while its native display is broken.
  for name in ("refinement_validation_spheres",):
    if defaults.get(name) is not None:
      _set_refinement_validation_display(name, defaults[name], show_dialog=False)
  return False


def _schedule_persistent_display_defaults_reapply():
  """Reapply saved defaults after Coot's own startup/state restoration."""
  glib_module = _coot_gui_repository_module("GLib")
  if glib_module is None:
    print("WARNING:: GLib unavailable; delayed coot_trimmings defaults were not scheduled")
    return False
  try:
    # Reassert at more than one point because state restoration timing varies
    # between Coot builds.  Each callback runs once (it returns False).
    source_ids = [
      glib_module.timeout_add(250, _apply_persistent_display_defaults),
      glib_module.timeout_add(1000, _apply_persistent_display_defaults),
      glib_module.timeout_add(3000, _apply_persistent_display_defaults),
    ]
    print("coot_trimmings: scheduled deferred defaults via GLib:", source_ids)
    return True
  except Exception as error:
    print(f"WARNING:: could not schedule display-default restoration: {error}")
    traceback.print_exc()
    return False


_apply_persistent_display_defaults()
_schedule_persistent_display_defaults_reapply()

def jiggle_fit_active_non_polymer_residue():
  residue = _active_residue_or_status()
  if not residue:
    return None
  mol_id = residue[0]
  ch_id = residue[1]
  resno = residue[2]
  ins_code = residue[3]
  if _residue_is_polymer(mol_id, ch_id, resno, ins_code):
    add_status_bar_text("Jiggle fit only applies to non-polymer residues")
    return None
  add_status_bar_text("Jiggle fitting active non-polymer residue")
  return fit_to_map_by_random_jiggle(mol_id, ch_id, resno, ins_code, 100, 1.0)


def smart_copy_active_non_polymer_residue():
  global SMART_COPY_SOURCE_CENTRE
  global SMART_COPY_RESIDUE_NAME
  residue = _active_residue_or_status()
  if not residue:
    return None
  mol_id = residue[0]
  ch_id = residue[1]
  resno = residue[2]
  ins_code = residue[3]
  if _residue_is_polymer(mol_id, ch_id, resno, ins_code):
    add_status_bar_text("Active residue is polymer; nothing copied")
    return None
  atom_selection = _smart_copy_atom_selection(ch_id, resno, ins_code)
  if atom_selection is None:
    add_status_bar_text("Smart copy does not yet support insertion-code residues")
    return None
  view_state = _capture_view_state()
  recentre_state = recentre_on_read_pdb()
  copied_imol = -1
  try:
    set_recentre_on_read_pdb(0)
    source_centre = residue_centre_py(mol_id, ch_id, resno, ins_code)
    source_name = residue_name(mol_id, ch_id, resno, ins_code)
    if not source_name or not isinstance(source_centre, (list, tuple)) or len(source_centre) != 3:
      raise ValueError("Copied residue information is unavailable")
    source_centre = tuple(float(value) for value in source_centre)
    if not all(math.isfinite(value) for value in source_centre):
      raise ValueError("Copied residue centre is invalid")
    copied_imol = new_molecule_by_atom_selection(mol_id, atom_selection)
    if copied_imol == -1 or coot.n_residues(copied_imol) != 1:
      raise ValueError("Unable to copy active residue")
  except Exception as error:
    if copied_imol != -1 and copied_imol in molecule_number_list():
      close_molecule(copied_imol)
    add_status_bar_text("Smart copy failed; previous copy retained: {0}".format(error))
    return None
  finally:
    set_recentre_on_read_pdb(recentre_state)
    _restore_view_state(view_state)
  _close_smart_copy_template()
  _set_smart_copy_template(copied_imol)
  SMART_COPY_SOURCE_CENTRE = source_centre
  SMART_COPY_RESIDUE_NAME = source_name
  add_status_bar_text("Copied {resname} for smart paste".format(resname=SMART_COPY_RESIDUE_NAME))


def _smart_paste_target_molecule():
  residue = active_residue()
  if residue:
    return residue[0]
  target_mol_id = go_to_atom_molecule_number()
  if target_mol_id in model_molecule_list():
    return target_mol_id
  return None


def _smart_paste_target_label(target_mol_id):
  return _molecule_id_name_label(target_mol_id)


def _prepare_smart_paste_molecule():
  if SMART_COPY_TEMPLATE_IMOL not in molecule_number_list():
    add_status_bar_text("No copied ligand, ion, or water available")
    return None
  if SMART_COPY_SOURCE_CENTRE is None:
    add_status_bar_text("Copied residue centre is unavailable")
    return None
  paste_imol = copy_molecule(SMART_COPY_TEMPLATE_IMOL)
  if paste_imol == -1:
    add_status_bar_text("Unable to prepare pasted residue")
    return None
  pointer_position = _rotation_centre_xyz()
  dx = pointer_position[0] - SMART_COPY_SOURCE_CENTRE[0]
  dy = pointer_position[1] - SMART_COPY_SOURCE_CENTRE[1]
  dz = pointer_position[2] - SMART_COPY_SOURCE_CENTRE[2]
  translate_molecule_by(paste_imol, dx, dy, dz)
  return paste_imol


def _smart_paste_new_molecule():
  paste_imol = _prepare_smart_paste_molecule()
  if paste_imol is None:
    return None
  set_mol_displayed(paste_imol, 1)
  set_mol_active(paste_imol, 1)
  add_status_bar_text(
    "Pasted {resname} as a new molecule ({target}); ready to paste again".format(
      resname=(SMART_COPY_RESIDUE_NAME or "copied residue"),
      target=_smart_paste_target_label(paste_imol),
    )
  )
  return paste_imol


def _smart_paste_into_molecule(target_mol_id):
  paste_imol = _prepare_smart_paste_molecule()
  if paste_imol is None:
    return None
  try:
    before = coot.n_residues(target_mol_id)
    added = coot.n_residues(paste_imol)
    if before < 0 or added != 1:
      add_status_bar_text("Unable to merge pasted residue; copy retained")
      return None
    # Coot's homogeneous-residue merge can add a residue but return [0].
    # n_atoms() excludes HETATM records, so verify the residue count instead.
    merge_molecules_py([paste_imol], target_mol_id)
    if coot.n_residues(target_mol_id) != before + added:
      add_status_bar_text("Unable to confirm pasted residue in active model; copy retained")
      return None
  finally:
    if paste_imol in molecule_number_list():
      close_molecule(paste_imol)
  _mark_model_changed(target_mol_id)
  add_status_bar_text(
    "Pasted {resname} at pointer in {target}; ready to paste again".format(
      resname=(SMART_COPY_RESIDUE_NAME or "copied residue"),
      target=_smart_paste_target_label(target_mol_id),
    )
  )
  return target_mol_id


def smart_paste_copied_non_polymer_residue():
  if SMART_COPY_TEMPLATE_IMOL not in molecule_number_list():
    add_status_bar_text("No copied ligand, ion, or water available")
    return None
  target_mol_id = _smart_paste_target_molecule()
  if target_mol_id is None:
    return _smart_paste_new_molecule()
  if "Gtk" not in globals():
    return _smart_paste_into_molecule(target_mol_id)
  generic_confirm_dialog(
    "Smart paste copied residue",
    "Do you want to paste into the active molecule "
    f"({_smart_paste_target_label(target_mol_id)}), or create a new molecule?",
    "Create a new molecule",
    lambda: _smart_paste_new_molecule(),
    "Paste into active",
    lambda: _smart_paste_into_molecule(target_mol_id),
  )


def go_to_nearest_density_peak():
  map_id = _scrollable_map_or_status()
  if map_id is None:
    return None
  centre = _rotation_centre_xyz()
  def density_here(point):
    try:
      return density_at_point(map_id, point[0], point[1], point[2])
    except Exception:
      return None

  peak_point = [centre[0], centre[1], centre[2]]
  peak_density = density_here(peak_point)
  if peak_density is None:
    add_status_bar_text("Unable to evaluate map density at the current position")
    return None

  for step_size in [0.5, 0.25, 0.1, 0.05]:
    for _step_index in range(40):
      best_point = peak_point
      best_density = peak_density
      for dx in [-step_size, 0.0, step_size]:
        for dy in [-step_size, 0.0, step_size]:
          for dz in [-step_size, 0.0, step_size]:
            if dx == 0.0 and dy == 0.0 and dz == 0.0:
              continue
            trial_point = [peak_point[0] + dx, peak_point[1] + dy, peak_point[2] + dz]
            trial_density = density_here(trial_point)
            if trial_density is None:
              continue
            if trial_density > best_density:
              best_point = trial_point
              best_density = trial_density
      if best_point == peak_point:
        break
      peak_point = best_point
      peak_density = best_density

  if _distance_sq(peak_point, centre) > 4.0 * 4.0:
    add_status_bar_text("No nearby density peak within 4 A")
    return None

  set_rotation_centre(peak_point[0], peak_point[1], peak_point[2])
  return peak_point


def _active_analysis_map_or_status():
  """Return the currently active map for density-analysis helpers."""
  map_id = scroll_wheel_map()
  if map_id != -1 and map_id in map_molecule_list():
    return map_id
  map_id = imol_refinement_map()
  if map_id != -1 and map_id in map_molecule_list():
    return map_id
  add_status_bar_text("No active map")
  return None


def _active_polymer_molecule_or_status():
  """Resolve the current model molecule for polymer-scoped validation helpers."""
  reference = _navigation_reference_residue()
  if reference is None:
    return None
  return reference["mol_id"]


def _residue_atom_records_and_xyz(imol, chain_id, resno, ins_code, snapshot=None):
  """Parse a residue once into both atom records and an atom-name lookup map."""
  atom_info = (snapshot.get((chain_id, resno, ins_code or ""), {}).get("atoms", [])
               if snapshot is not None else residue_info_py(imol, chain_id, resno, ins_code))
  if not isinstance(atom_info, list):
    return [], {}

  atom_records = []
  atom_xyz = {}
  for atom in atom_info:
    try:
      atom_name = atom[0][0].strip()
      element = atom[1][2].strip().upper()
      xyz = list(atom[2])
    except Exception:
      continue
    if atom_name and len(xyz) == 3:
      atom_records.append({
        "name": atom_name,
        "element": element,
        "xyz": xyz,
      })
      atom_xyz.setdefault(atom_name, xyz)
  return atom_records, atom_xyz


def _residue_atom_records(imol, chain_id, resno, ins_code):
  """Return parsed atom records with names, coordinates, and elements."""
  atom_records, _atom_xyz = _residue_atom_records_and_xyz(imol, chain_id, resno, ins_code)
  return atom_records


def _pseudo_cb_direction_from_backbone(atom_xyz):
  """Construct an alanine-like pseudo-CB direction from backbone geometry."""
  if not all(atom_name in atom_xyz for atom_name in ("N", "CA", "C")):
    return None
  anchor = atom_xyz["CA"]
  n_vector = _normalize_vector(_vector_subtract(atom_xyz["N"], anchor))
  c_vector = _normalize_vector(_vector_subtract(atom_xyz["C"], anchor))
  if n_vector is None or c_vector is None:
    return None

  hydrogen_midpoint = [
    -(n_vector[0] + c_vector[0]) * 0.5,
    -(n_vector[1] + c_vector[1]) * 0.5,
    -(n_vector[2] + c_vector[2]) * 0.5,
  ]
  plane_normal = _normalize_vector(_vector_cross(n_vector, c_vector))
  if plane_normal is None:
    return None
  if "O" in atom_xyz:
    o_vector = _normalize_vector(_vector_subtract(atom_xyz["O"], anchor))
    if o_vector is not None and _vector_dot(plane_normal, o_vector) > 0.0:
      plane_normal = _scale_vector(plane_normal, -1.0)

  perpendicular_scale_sq = max(0.0, 1.0 - _vector_dot(hydrogen_midpoint, hydrogen_midpoint))
  perpendicular_scale = math.sqrt(perpendicular_scale_sq)
  return _normalize_vector([
    hydrogen_midpoint[0] + perpendicular_scale * plane_normal[0],
    hydrogen_midpoint[1] + perpendicular_scale * plane_normal[1],
    hydrogen_midpoint[2] + perpendicular_scale * plane_normal[2],
  ])


def _emringer_ring_frame(atom_xyz):
  """Build the CA-CB frame used to sweep a virtual CG around chi1."""
  if not all(atom_name in atom_xyz for atom_name in ("N", "CA", "C")):
    return None

  ca_xyz = atom_xyz["CA"]
  if "CB" in atom_xyz:
    cb_xyz = atom_xyz["CB"]
    axis_unit = _normalize_vector(_vector_subtract(cb_xyz, ca_xyz))
  else:
    pseudo_cb_direction = _pseudo_cb_direction_from_backbone(atom_xyz)
    if pseudo_cb_direction is None:
      return None
    axis_unit = pseudo_cb_direction
    cb_xyz = [
      ca_xyz[0] + EMRINGER_HELPER_CA_CB_BOND_LENGTH * axis_unit[0],
      ca_xyz[1] + EMRINGER_HELPER_CA_CB_BOND_LENGTH * axis_unit[1],
      ca_xyz[2] + EMRINGER_HELPER_CA_CB_BOND_LENGTH * axis_unit[2],
    ]

  if axis_unit is None:
    return None

  ring_u = _project_vector_perpendicular(_vector_subtract(atom_xyz["N"], ca_xyz), axis_unit)
  if ring_u is None:
    ring_u = _project_vector_perpendicular(_vector_subtract(atom_xyz["C"], ca_xyz), axis_unit)
  if ring_u is None:
    return None
  ring_v = _normalize_vector(_vector_cross(axis_unit, ring_u))
  if ring_v is None:
    return None

  return {
    "ca_xyz": ca_xyz,
    "cb_xyz": cb_xyz,
    "anchor_xyz": cb_xyz,
    "axis_unit": axis_unit,
    "ring_u": ring_u,
    "ring_v": ring_v,
  }


def _emringer_torsion_ring_frame(atom_xyz, reference_atom_name, axis_start_atom_name, axis_end_atom_name):
  """Build a generic torsion frame for sweeping a downstream sidechain atom."""
  if not all(atom_name in atom_xyz for atom_name in (reference_atom_name, axis_start_atom_name, axis_end_atom_name)):
    return None

  axis_start_xyz = atom_xyz[axis_start_atom_name]
  axis_end_xyz = atom_xyz[axis_end_atom_name]
  axis_unit = _normalize_vector(_vector_subtract(axis_end_xyz, axis_start_xyz))
  if axis_unit is None:
    return None

  ring_u = _project_vector_perpendicular(
    _vector_subtract(atom_xyz[reference_atom_name], axis_start_xyz),
    axis_unit,
  )
  if ring_u is None:
    return None
  ring_v = _normalize_vector(_vector_cross(axis_unit, ring_u))
  if ring_v is None:
    return None

  return {
    "anchor_xyz": axis_end_xyz,
    "axis_unit": axis_unit,
    "ring_u": ring_u,
    "ring_v": ring_v,
  }


def _virtual_torsion_point(ring_frame, bond_length, angle_degrees):
  anchor_xyz = ring_frame["anchor_xyz"]
  axis_unit = ring_frame["axis_unit"]
  ring_u = ring_frame["ring_u"]
  ring_v = ring_frame["ring_v"]
  axial_offset = bond_length / 3.0
  radial_offset = bond_length * math.sqrt(8.0 / 9.0)
  angle_radians = math.radians(angle_degrees)
  ring_direction = [
    math.cos(angle_radians) * ring_u[0] + math.sin(angle_radians) * ring_v[0],
    math.cos(angle_radians) * ring_u[1] + math.sin(angle_radians) * ring_v[1],
    math.cos(angle_radians) * ring_u[2] + math.sin(angle_radians) * ring_v[2],
  ]
  return [
    anchor_xyz[0] + axial_offset * axis_unit[0] + radial_offset * ring_direction[0],
    anchor_xyz[1] + axial_offset * axis_unit[1] + radial_offset * ring_direction[1],
    anchor_xyz[2] + axial_offset * axis_unit[2] + radial_offset * ring_direction[2],
  ]


def _sample_virtual_stage_density_peak(map_id, ring_frame, bond_length,
                                       step_degrees=EMRINGER_HELPER_STEP_DEGREES,
                                       fine_step_degrees=EMRINGER_HELPER_FINE_STEP_DEGREES):
  """Find the strongest density peak along a virtual sidechain torsion sweep."""
  if step_degrees <= 0.0:
    step_degrees = EMRINGER_HELPER_STEP_DEGREES

  best_density = None
  best_angle = None
  best_point = None
  coarse_density_sum = 0.0
  coarse_density_count = 0
  n_steps = max(1, int(round(360.0 / step_degrees)))
  for step_index in range(n_steps):
    angle_degrees = (360.0 * step_index) / n_steps
    sample_point = _virtual_torsion_point(ring_frame, bond_length, angle_degrees)
    density_value = density_at_point(map_id, sample_point[0], sample_point[1], sample_point[2])
    coarse_density_sum += density_value
    coarse_density_count += 1
    if best_density is None or density_value > best_density:
      best_density = density_value
      best_angle = angle_degrees
      best_point = sample_point

  if best_density is None:
    return None
  ring_mean_density = coarse_density_sum / float(coarse_density_count)

  if fine_step_degrees > 0.0 and fine_step_degrees < step_degrees and best_angle is not None:
    search_start = best_angle - step_degrees
    search_stop = best_angle + step_degrees
    n_fine_steps = max(1, int(round((search_stop - search_start) / fine_step_degrees)))
    for step_index in range(n_fine_steps + 1):
      angle_degrees = search_start + step_index * fine_step_degrees
      sample_point = _virtual_torsion_point(ring_frame, bond_length, angle_degrees)
      density_value = density_at_point(map_id, sample_point[0], sample_point[1], sample_point[2])
      if density_value > best_density:
        best_density = density_value
        best_angle = angle_degrees
        best_point = sample_point

  return {
    "density": best_density,
    "angle_degrees": best_angle % 360.0 if best_angle is not None else None,
    "point": best_point,
    "ring_mean_density": ring_mean_density,
  }


def _sample_virtual_cg_density_peak(map_id, ring_frame, step_degrees=EMRINGER_HELPER_STEP_DEGREES):
  """Find the strongest density peak along the virtual-CG chi1 sweep."""
  return _sample_virtual_stage_density_peak(
    map_id,
    ring_frame,
    EMRINGER_HELPER_CB_CG_BOND_LENGTH,
    step_degrees=step_degrees,
  )


def _emringer_late_stage_definitions(residue_name_here):
  """Return later-torsion sweeps for longer, flexible sidechains."""
  return EMRINGER_HELPER_LATE_STAGE_DEFINITIONS.get(residue_name_here, [])


def _emringer_chi1_blocker_atom_names(residue_name_here):
  """Return the atom names that genuinely occupy the chi1/gamma position.

  The first-pass helper used every heavy atom beyond CB as an occupancy check,
  which was too blunt: a distal atom could hide a real gamma-position outlier.
  For the chi1 sweep we only want atoms that belong at the gamma stage for the
  residue, plus a conservative fallback for unfamiliar components.
  """
  return EMRINGER_HELPER_CHI1_BLOCKER_ATOM_NAMES.get(residue_name_here)


def _emringer_chi1_blocker_atoms(atom_records, residue_name_here):
  blocker_atom_names = _emringer_chi1_blocker_atom_names(residue_name_here)
  blocker_atoms = []
  for atom in atom_records:
    if atom["element"] in ("H", "D"):
      continue
    atom_name = atom["name"]
    if blocker_atom_names is None:
      if atom_name in EMRINGER_HELPER_BACKBONE_ATOM_NAMES:
        continue
    else:
      if atom_name not in blocker_atom_names:
        continue
    blocker_atoms.append(atom)
  return blocker_atoms


def _emringer_blocker_atoms(atom_records, blocker_atom_names):
  blocker_atoms = []
  for atom in atom_records:
    if atom["element"] in ("H", "D"):
      continue
    if atom["name"] in blocker_atom_names:
      blocker_atoms.append(atom)
  return blocker_atoms


def _nucleotide_base_atoms(atom_records, glycosidic_atom_name=None):
  """Return heavy atoms that belong to the base rather than the sugar/phosphate."""
  base_atoms = []
  for atom in atom_records:
    atom_name = atom["name"]
    if atom["element"] in ("H", "D"):
      continue
    if "'" in atom_name or atom_name.startswith("OP") or atom_name in {"P", "O1P", "O2P", "O3P"}:
      continue
    if glycosidic_atom_name and atom_name == glycosidic_atom_name:
      continue
    base_atoms.append(atom)
  return base_atoms


def _nucleotide_stage_specs(atom_xyz, atom_records):
  """Build a glycosidic torsion sweep for nucleic acid bases."""
  if "C1'" not in atom_xyz:
    return []

  reference_atom_name = next((name for name in ("O4'", "C2'", "C3'") if name in atom_xyz), None)
  if reference_atom_name is None:
    return []

  if "N9" in atom_xyz:
    glycosidic_atom_name = "N9"
    blocker_atom_names = {"C4", "C8"}
  elif "N1" in atom_xyz:
    glycosidic_atom_name = "N1"
    blocker_atom_names = {"C2", "C6"}
  else:
    return []

  ring_frame = _emringer_torsion_ring_frame(atom_xyz, reference_atom_name, "C1'", glycosidic_atom_name)
  if ring_frame is None:
    return []

  blocker_atoms = _emringer_blocker_atoms(atom_records, blocker_atom_names)
  if not blocker_atoms:
    blocker_atoms = _nucleotide_base_atoms(atom_records, glycosidic_atom_name)
  if not blocker_atoms:
    return []

  return [{
    "label": "chi",
    "ring_frame": ring_frame,
    "bond_length": 1.38,
    "blocker_atoms": blocker_atoms,
  }]


def _emringer_stage_specs(atom_xyz, atom_records, residue_name_here):
  """Build the torsion stages that should be tested for a residue."""
  if "C1'" in atom_xyz:
    return _nucleotide_stage_specs(atom_xyz, atom_records)

  stage_specs = []

  chi1_ring_frame = _emringer_ring_frame(atom_xyz)
  if chi1_ring_frame is not None:
    stage_specs.append({
      "label": "chi1",
      "ring_frame": chi1_ring_frame,
      "bond_length": EMRINGER_HELPER_CB_CG_BOND_LENGTH,
      "blocker_atoms": _emringer_chi1_blocker_atoms(atom_records, residue_name_here),
    })

  for stage_definition in _emringer_late_stage_definitions(residue_name_here):
    ring_frame = _emringer_torsion_ring_frame(
      atom_xyz,
      stage_definition["reference"],
      stage_definition["axis_start"],
      stage_definition["axis_end"],
    )
    if ring_frame is None:
      continue
    stage_specs.append({
      "label": stage_definition["label"],
      "ring_frame": ring_frame,
      "bond_length": stage_definition["bond_length"],
      "blocker_atoms": _emringer_blocker_atoms(atom_records, stage_definition["blockers"]),
    })

  return stage_specs


def _nearest_stage_blocker(stage_spec, peak_point):
  """Return the nearest blocker atom to a stage peak, if any."""
  nearest_distance_sq = None
  nearest_blocker_name = None
  for atom in stage_spec["blocker_atoms"]:
    distance_sq = _distance_sq(atom["xyz"], peak_point)
    if nearest_distance_sq is None or distance_sq < nearest_distance_sq:
      nearest_distance_sq = distance_sq
      nearest_blocker_name = atom["name"]
  return nearest_blocker_name, nearest_distance_sq


def _mean_density_at_atoms(map_id, atom_records):
  """Average map density at the current placed stage atoms, if any exist."""
  if not atom_records:
    return None

  density_sum = 0.0
  density_count = 0
  for atom in atom_records:
    xyz = atom["xyz"]
    density_sum += density_at_point(map_id, xyz[0], xyz[1], xyz[2])
    density_count += 1
  if density_count == 0:
    return None
  return density_sum / float(density_count)


def _backbone_support_positions(atom_xyz):
  """Return the local backbone/sugar positions used to prefilter weak residues."""
  if "C1'" in atom_xyz or "P" in atom_xyz:
    atom_names = ("P", "O5'", "C5'", "C4'", "C3'", "O3'", "C1'")
  else:
    atom_names = ("N", "CA", "C", "O")
  return [atom_xyz[atom_name] for atom_name in atom_names if atom_name in atom_xyz]


def _backbone_density_support(map_id, atom_xyz, density_threshold):
  """Return the fraction of backbone support points that are above threshold."""
  backbone_positions = _backbone_support_positions(atom_xyz)
  if not backbone_positions:
    return 0.0, 0, 0

  supported_points = 0
  for position in backbone_positions:
    if density_at_point(map_id, position[0], position[1], position[2]) >= density_threshold:
      supported_points += 1
  total_points = len(backbone_positions)
  return (float(supported_points) / total_points), supported_points, total_points


def _residue_display_point(atom_records, atom_xyz):
  """Choose a useful jump target for non-polymer residues and summaries."""
  if atom_records:
    xs = [atom["xyz"][0] for atom in atom_records]
    ys = [atom["xyz"][1] for atom in atom_records]
    zs = [atom["xyz"][2] for atom in atom_records]
    return [sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs)]
  if atom_xyz:
    return next(iter(atom_xyz.values()))
  return None


def _evaluate_emringer_stage(map_id, stage_spec):
  """Evaluate one torsion stage and return its peak/current-density summary."""
  peak = _sample_virtual_stage_density_peak(map_id, stage_spec["ring_frame"], stage_spec["bond_length"])
  if peak is None:
    return None

  peak_angle = peak["angle_degrees"]
  max_local_support_density = None
  support_density_sum = peak["density"]
  support_density_count = 1
  if peak_angle is not None and EMRINGER_HELPER_MISFIT_SUPPORT_STEP_DEGREES > 0.0:
    for angle_offset in (-EMRINGER_HELPER_MISFIT_SUPPORT_STEP_DEGREES, EMRINGER_HELPER_MISFIT_SUPPORT_STEP_DEGREES):
      support_point = _virtual_torsion_point(
        stage_spec["ring_frame"],
        stage_spec["bond_length"],
        peak_angle + angle_offset,
      )
      support_density = density_at_point(map_id, support_point[0], support_point[1], support_point[2])
      if max_local_support_density is None or support_density > max_local_support_density:
        max_local_support_density = support_density
      support_density_sum += support_density
      support_density_count += 1

  current_stage_density = _mean_density_at_atoms(map_id, stage_spec["blocker_atoms"])
  nearest_blocker_name, nearest_distance_sq = _nearest_stage_blocker(stage_spec, peak["point"])
  density_gain = peak["density"] - current_stage_density if current_stage_density is not None else peak["density"]
  ring_prominence = peak["density"] - peak["ring_mean_density"]
  mean_triplet_density = support_density_sum / float(support_density_count)
  return {
    "label": stage_spec["label"],
    "bond_length": stage_spec["bond_length"],
    "ring_frame": stage_spec["ring_frame"],
    "peak": peak,
    "current_stage_density": current_stage_density,
    "nearest_blocker_name": nearest_blocker_name,
    "nearest_distance_sq": nearest_distance_sq,
    "density_gain": density_gain,
    "ring_prominence": ring_prominence,
    "max_local_support_density": max_local_support_density,
    "mean_triplet_density": mean_triplet_density,
  }


def _weak_stage_inset_peak_density(map_id, stage_result, inset_buffer=EMRINGER_HELPER_WEAK_STAGE_INSET_BUFFER):
  """Sample a slightly inset torsion ring to avoid flagging borderline weak cases."""
  inset_bond_length = stage_result["bond_length"] - inset_buffer
  if inset_bond_length <= 0.0:
    return stage_result["peak"]["density"]
  inset_peak = _sample_virtual_stage_density_peak(
    map_id,
    stage_result["ring_frame"],
    inset_bond_length,
  )
  if inset_peak is None:
    return stage_result["peak"]["density"]
  return max(stage_result["peak"]["density"], inset_peak["density"])


def _heavy_atom_density_summary(map_id, atom_records, density_threshold=None):
  """Sample heavy-atom density once and return the summary statistics."""
  heavy_atoms = []
  max_density = None
  max_atom = None
  density_sum = 0.0
  outside_count = 0

  for atom in atom_records:
    if atom["element"] in ("H", "D"):
      continue
    heavy_atoms.append(atom)
    x, y, z = atom["xyz"]
    density_value = density_at_point(map_id, x, y, z)
    density_sum += density_value
    if max_density is None or density_value > max_density:
      max_density = density_value
      max_atom = atom
    if density_threshold is not None and density_value < density_threshold:
      outside_count += 1

  heavy_atom_count = len(heavy_atoms)
  mean_density = density_sum / float(heavy_atom_count) if heavy_atom_count else None
  return {
    "heavy_atoms": heavy_atoms,
    "count": heavy_atom_count,
    "max_density": max_density,
    "max_atom": max_atom,
    "mean_density": mean_density,
    "outside_count": outside_count,
  }


def _format_residue_id(chain_id, resno, ins_code):
  return f"{chain_id}{resno}{ins_code or ''}"


def _display_residue_name(residue_name):
  """Return a readable residue label for GUI lists."""
  if not residue_name:
    return ""
  if residue_name in POLYMER_RESIDUE_NAMES and len(residue_name) == 3 and residue_name.isalpha():
    return residue_name.capitalize()
  return residue_name.upper()


def _odd_residue_dialog_label(chain_id, resno, ins_code, residue_name=None):
  """Compact GUI label: chain:resno[icode] residue-type."""
  residue_id = f"{chain_id}:{resno}{ins_code or ''}"
  if residue_name:
    return f"{residue_id} {_display_residue_name(residue_name)}"
  return residue_id


def _log_odd_residue_results(mol_id, map_id, peak_threshold, categorized_details):
  """Print full odd-residue diagnostics to the log in dialog order."""
  print(
    "Odd residues for molecule #{0} at the current threshold of the scrollable map "
    "(Map #{1}); threshold={2:.4f}".format(mol_id, map_id, peak_threshold)
  )
  for category_name, detail_lines in categorized_details:
    print(f"{category_name} ({len(detail_lines)})")
    for detail_line in detail_lines:
      print(f"  {detail_line}")


def _append_odd_residue_entry(categorized_entries, category_name, sort_key, dialog_label, point, detail_label,
                              navigation_metadata=None):
  """Store one GUI/log entry tuple in the requested odd-residue category."""
  gui_entry = [dialog_label, point[0], point[1], point[2]]
  if navigation_metadata is not None:
    gui_entry.append(navigation_metadata)
  categorized_entries[category_name].append((sort_key, gui_entry, detail_label))


def _sorted_categorized_outputs(categorized_entries, category_order=None):
  """Sort each category once, then split into GUI entries and log-detail lines."""
  gui_categories = []
  detail_categories = []
  ordered_categories = category_order or tuple(categorized_entries.keys())
  for category_name in ordered_categories:
    sorted_entries = sorted(categorized_entries[category_name], key=lambda item: item[0])
    gui_categories.append((category_name, [entry for _sort_key, entry, _detail in sorted_entries]))
    detail_categories.append((category_name, [detail for _sort_key, _entry, detail in sorted_entries]))
  return gui_categories, detail_categories


def _odd_residue_filter_categories(gui_categories, chain_id=None):
  filtered = [(category, [entry for entry in entries
              if chain_id is None or (len(entry) > 4 and isinstance(entry[4], dict)
                                      and entry[4].get("chain_id") == chain_id)])
              for category, entries in gui_categories]
  return _numbered_gui_categories(filtered)


def _numbered_gui_categories(gui_categories):
  """Prefix categorized hit labels with per-category numbers for easier navigation."""
  numbered_categories = []
  for category_name, entries in gui_categories:
    numbered_entries = []
    for entry_index, entry in enumerate(entries, start=1):
      if not isinstance(entry, list) or not entry:
        numbered_entries.append(entry)
        continue
      numbered_entry = list(entry)
      numbered_entry[0] = f"{entry_index}. {numbered_entry[0]}"
      numbered_entries.append(numbered_entry)
    numbered_categories.append((category_name, numbered_entries))
  return numbered_categories


def _odd_residue_navigation_metadata(mol_id, chain_id, resno, ins_code, serial_number):
  """Describe a polymer residue jump target for the smart navigation helper."""
  return {
    "type": "polymer_residue",
    "mol_id": mol_id,
    "chain_id": chain_id,
    "resno": resno,
    "ins_code": ins_code or "",
    "serial_number": serial_number,
  }


def _odd_residue_status_suffix(category_name):
  """Return a concise category label for status-bar navigation messages."""
  status_suffixes = {
    "Possible Misfits": "Possible Misfit?",
    "Weak Sidechains": "Weak Sidechain",
    "Weak Backbone": "Weak Backbone",
    "Weak Waters/Ions": "Weak Water/Ion",
    "Weak Ligands": "Weak Ligand",
  }
  return status_suffixes.get(category_name, category_name)


def _missing_atom_residue_keys(mol_id):
  """Return residue keys flagged by Coot as having missing atoms."""
  try:
    missing_specs = missing_atom_info(mol_id)
  except Exception:
    try:
      missing_specs = missing_atom_info_py(mol_id)
    except Exception:
      return set()

  residue_keys = set()
  if not isinstance(missing_specs, (list, tuple)):
    return residue_keys

  for spec in missing_specs:
    if not isinstance(spec, (list, tuple)) or len(spec) < 2:
      continue
    chain_id = spec[0]
    resno = spec[1]
    ins_code = spec[2] if len(spec) > 2 else ""
    if isinstance(chain_id, str):
      residue_keys.add((chain_id, resno, ins_code or ""))
  return residue_keys


def _odd_residue_peak_threshold_or_status(map_id):
  """Return the current positive contour level used as the odd-residues threshold."""
  contour_level = get_contour_level_absolute(map_id)
  if contour_level <= 0.0:
    info_dialog(
      "Odd residues requires a positive scrollable-map contour level.\n"
      f"Coot reported {contour_level:.4f} for map #{map_id}."
    )
    return None
  return contour_level


def _odd_residue_dialog_title(mol_id, map_id):
  """Return the odd-residues dialog title with current map level details."""
  contour_level = get_contour_level_absolute(map_id)
  try:
    contour_sigma = get_contour_level_in_sigma(map_id)
  except Exception:
    contour_sigma = None

  if contour_sigma is not None and contour_sigma is not False:
    return (
      "Odd residues in molecule #{0} at the current threshold of the scrollable map\n"
      "(Map #{1} at level {2:.3f}, {3:.1f} sigma):"
    ).format(mol_id, map_id, contour_level, contour_sigma)

  return (
    "Odd residues in molecule #{0} at the current threshold of the scrollable map\n"
    "(Map #{1} at level {2:.3f}):"
  ).format(mol_id, map_id, contour_level)


def _stage_result_has_local_support(stage_result, peak_threshold):
  """Require the misfit peak to be supported by nearby torsion samples."""
  max_local_support_density = stage_result["max_local_support_density"]
  return (
    (max_local_support_density is not None and max_local_support_density >= peak_threshold)
    or stage_result["mean_triplet_density"] >= peak_threshold
  )


def _stage_result_is_possible_misfit(stage_result, peak_threshold, misfit_peak_threshold,
                                     misfit_gain_threshold, occupied_distance_sq):
  """Return True if this torsion-stage peak looks like a genuine possible misfit."""
  current_stage_density = stage_result["current_stage_density"]
  nearest_distance_sq = stage_result["nearest_distance_sq"]
  prominence_threshold = peak_threshold * EMRINGER_HELPER_MISFIT_MIN_PROMINENCE_TO_CONTOUR_RATIO
  if stage_result["peak"]["density"] < misfit_peak_threshold:
    return False
  if stage_result["ring_prominence"] < prominence_threshold:
    return False
  if not _stage_result_has_local_support(stage_result, peak_threshold):
    return False
  if current_stage_density is not None and stage_result["peak"]["density"] <= current_stage_density:
    return False
  if current_stage_density is not None and stage_result["density_gain"] < misfit_gain_threshold:
    return False
  if nearest_distance_sq is not None and nearest_distance_sq <= occupied_distance_sq:
    return False
  return True


def _possible_misfit_detail_label(residue_id, residue_name_here, stage_result, peak_threshold):
  """Format the log detail for one possible-misfit stage hit."""
  contour_ratio = stage_result["peak"]["density"] / peak_threshold if peak_threshold > 0.0 else 0.0
  if stage_result["nearest_distance_sq"] is None:
    placement_note = f"no {stage_result['label']}-stage atoms placed"
  else:
    placement_note = "{0} is {1:.1f} A away".format(
      stage_result["nearest_blocker_name"] or f"nearest {stage_result['label']}-stage atom",
      math.sqrt(stage_result["nearest_distance_sq"]),
    )
  if stage_result["current_stage_density"] is None:
    density_note = "current stage unbuilt"
  else:
    density_note = "current stage density {0:.4f}, gain {1:.4f}".format(
      stage_result["current_stage_density"],
      stage_result["density_gain"],
    )
  return (
    f"{residue_id} {residue_name_here}: peak {stage_result['peak']['density']:.4f} "
    f"({contour_ratio:.1f}x threshold) at {stage_result['label']} {stage_result['peak']['angle_degrees']:.0f} deg; "
    f"{density_note}; ring mean {stage_result['peak']['ring_mean_density']:.4f}, "
    f"prominence {stage_result['ring_prominence']:.4f}; local support max {stage_result['max_local_support_density'] or 0.0:.4f}, "
    f"triplet mean {stage_result['mean_triplet_density']:.4f}; {placement_note}"
  )


def _scan_polymer_stages(map_id, residue_id, residue_name_here, stage_specs, peak_threshold,
                         misfit_peak_threshold, misfit_gain_threshold, occupied_distance_sq):
  """Evaluate all torsion stages for one polymer residue."""
  best_stage_result = None
  strongest_current_density = None
  strongest_buffered_peak_density = None
  has_stage_atoms = False
  possible_misfits = []

  for stage_index, stage_spec in enumerate(stage_specs):
    stage_result = _evaluate_emringer_stage(map_id, stage_spec)
    if stage_result is None:
      continue

    if best_stage_result is None or stage_result["peak"]["density"] > best_stage_result["peak"]["density"]:
      best_stage_result = stage_result

    current_stage_density = stage_result["current_stage_density"]
    if current_stage_density is not None and (
      strongest_current_density is None or current_stage_density > strongest_current_density
    ):
      strongest_current_density = current_stage_density

    buffered_peak_density = _weak_stage_inset_peak_density(map_id, stage_result)
    if strongest_buffered_peak_density is None or buffered_peak_density > strongest_buffered_peak_density:
      strongest_buffered_peak_density = buffered_peak_density

    if stage_result["nearest_distance_sq"] is not None:
      has_stage_atoms = True

    if _stage_result_is_possible_misfit(
      stage_result,
      peak_threshold,
      misfit_peak_threshold,
      misfit_gain_threshold,
      occupied_distance_sq,
    ):
      possible_misfits.append((stage_index, stage_result))

  return {
    "best_stage_result": best_stage_result,
    "strongest_current_density": strongest_current_density,
    "strongest_buffered_peak_density": strongest_buffered_peak_density,
    "has_stage_atoms": has_stage_atoms,
    "possible_misfits": possible_misfits,
  }


def _odd_dihedral_degrees(a, b, c, d):
  axis = _normalize_vector(_vector_subtract(c, b))
  if axis is None:
    return None
  v = _project_vector_perpendicular(_vector_subtract(a, b), axis)
  w = _project_vector_perpendicular(_vector_subtract(d, c), axis)
  if v is None or w is None:
    return None
  return math.degrees(math.atan2(_vector_dot(_vector_cross(axis, v), w), _vector_dot(v, w)))


def _odd_phi_psi_context_hits(angles, links):
  """Flag clear helix/strand switches inside a consistent two-residue context."""
  def region(pair):
    phi, psi = pair
    if -100 <= phi <= -30 and -80 <= psi <= -5:
      return "helix"
    if -180 <= phi <= -70 and (psi >= 70 or psi <= -150):
      return "strand"
    return None

  def clearly_opposite(pair, context):
    phi, psi = pair
    if context == "helix":
      return -180 <= phi <= -90 and (psi >= 90 or psi <= -160)
    return -95 <= phi <= -35 and -75 <= psi <= -5

  hits = []
  for index in range(3, len(angles)-3):
    local = angles[index-2:index+3]
    if any(pair is None or not all(math.isfinite(value) for value in pair) for pair in local):
      continue
    if not all(links[index-3:index+3]):
      continue
    neighbours = [region(local[i]) for i in (0, 1, 3, 4)]
    context = neighbours[0]
    if context is not None and all(value == context for value in neighbours) and clearly_opposite(local[2], context):
      hits.append((index, context, *local[2]))
  return hits


def _odd_backbone_geometry(records):
  usable = []
  for record in records:
    xyz = record["xyz"]
    usable.append(record["name"] in _REGISTER_PROTEIN_NAMES and
                  all(name in xyz for name in ("N", "CA", "C", "O")) and record["unambiguous"])
  links = [usable[i] and usable[i+1] and
           1.0**2 <= _distance_sq(records[i]["xyz"]["C"], records[i+1]["xyz"]["N"]) <= 1.8**2
           for i in range(len(records)-1)]
  angles = [None] * len(records)
  for i in range(1, len(records)-1):
    if not (links[i-1] and links[i]):
      continue
    xyz = records[i]["xyz"]
    phi = _odd_dihedral_degrees(records[i-1]["xyz"]["C"], xyz["N"], xyz["CA"], xyz["C"])
    psi = _odd_dihedral_degrees(xyz["N"], xyz["CA"], xyz["C"], records[i+1]["xyz"]["N"])
    if phi is not None and psi is not None:
      angles[i] = (phi, psi)
  return angles, links


def _odd_backbone_context_hits(records):
  angles, links = _odd_backbone_geometry(records)
  return _odd_phi_psi_context_hits(angles, links)


def _odd_backbone_hbond(records, acceptor, donor):
  """Approximate a peptide N-H direction from the preceding carbonyl."""
  if donor < 1 or records[donor]["name"] == "PRO":
    return False
  acceptor_xyz, donor_xyz, previous_xyz = (records[index]["xyz"] for index in (acceptor, donor, donor-1))
  direction = _normalize_vector(_vector_subtract(previous_xyz["C"], previous_xyz["O"]))
  if direction is None:
    return False
  hydrogen = tuple(donor_xyz["N"][axis] + direction[axis] for axis in range(3))
  def distance(a, b):
    return math.sqrt(_distance_sq(a, b))
  on, ch, oh, cn = (distance(a, b) for a, b in (
    (acceptor_xyz["O"], donor_xyz["N"]), (acceptor_xyz["C"], hydrogen),
    (acceptor_xyz["O"], hydrogen), (acceptor_xyz["C"], donor_xyz["N"])))
  if not all(math.isfinite(value) and value > 0.5 for value in (on, ch, oh, cn)) or not 2.3 <= on <= 3.5:
    return False
  energy = 27.888 * (1/on + 1/ch - 1/oh - 1/cn)
  return energy < -0.5 and oh < on


def _odd_ca_window_kind(ca, start):
  """Classify a C-alpha distance fingerprint over five consecutive residues."""
  def distance(offset):
    return math.sqrt(_distance_sq(ca[start], ca[start+offset]))
  d3, d4 = distance(3), distance(4)
  ca_torsion = _odd_dihedral_degrees(ca[start], ca[start+1], ca[start+2], ca[start+3])
  if ca_torsion is None:
    return None
  # RyR 3_10 segments can have a contracted i->i+3 distance while retaining
  # the characteristic i->i+4 span and torsion.
  if 5.35 <= d3 <= 6.35 and 7.5 <= d4 <= 9.2 and 63 <= ca_torsion <= 103:
    return "310"
  if 4.7 <= d3 <= 5.55 and 5.8 <= d4 <= 6.7 and 35 <= ca_torsion <= 65:
    return "alpha"
  if 5.5 <= d3 <= 6.2 and 4.6 <= d4 <= 5.2 and 15 <= ca_torsion <= 42:
    return "pi"
  return None


def _odd_helical_patterns(ca, links):
  """Infer tentative helical motifs from consecutive C-alpha distance fingerprints."""
  count = len(ca)
  kinds = [None] * max(0, count-4)
  for start in range(len(kinds)):
    if all(links[start:start+4]):
      kinds[start] = _odd_ca_window_kind(ca, start)
  features = []
  for motif, label in (("310", "Possible 3_10 helix"), ("pi", "Possible pi-bulge")):
    starts = [i for i, kind in enumerate(kinds)
              if kind == motif and (motif != "pi" or i+5 < count)]
    runs = []
    for start in starts:
      if runs and start == runs[-1][-1] + 1:
        runs[-1].append(start)
      else:
        runs.append([start])
    for run in runs:
      first, last = run[0], run[-1]
      flank_start = max(0, first-2)
      flank_end = min(len(kinds), last+3)
      left_flanks = kinds[flank_start:first]
      right_flanks = kinds[last+1:flank_end]
      alpha_flanks = (len(left_flanks) == 2 and len(right_flanks) == 2 and
                      all(kind == "alpha" for kind in left_flanks + right_flanks))
      if motif == "pi" and len(run) <= 2 and alpha_flanks:
        label_here = "Possible pi-bulge"
      elif motif == "310" and len(run) >= 2:
        label_here = label
      else:
        continue
      span = 3 if motif == "310" else 5
      features.append((first, last+max(4, span), label_here, len(run), span, tuple(run)))
  return sorted(features)


def _odd_helical_feature_hits(records):
  ca = [record["xyz"].get("CA") for record in records]
  usable = [record["name"] in _REGISTER_PROTEIN_NAMES and record["unambiguous"] and
            point is not None and len(point) == 3 and all(math.isfinite(v) for v in point)
            for record, point in zip(records, ca)]
  links = [usable[i] and usable[i+1] and 3.4**2 <= _distance_sq(ca[i], ca[i+1]) <= 4.1**2
           for i in range(len(records)-1)]
  safe_ca = [point if valid else (float("nan"),)*3 for point, valid in zip(ca, usable)]
  features = _odd_helical_patterns(safe_ca, links)
  supported_features = []
  for first, last, label, pattern_count, step, starts in features:
    hbond_flags = []
    for start in starts:
      acceptor, donor = start, start + step
      required = all(atom in records[acceptor]["xyz"] for atom in ("C", "O")) and all(
        atom in records[index]["xyz"] for index in (donor-1, donor) for atom in ("C", "O", "N"))
      peptide_path_ok = all(
        1.0**2 <= _distance_sq(records[index]["xyz"]["C"], records[index+1]["xyz"]["N"]) <= 1.8**2
        for index in range(acceptor, donor)
        if "C" in records[index]["xyz"] and "N" in records[index+1]["xyz"])
      peptide_path_ok = peptide_path_ok and all(
        "C" in records[index]["xyz"] and "N" in records[index+1]["xyz"]
        for index in range(acceptor, donor))
      if (not required or not peptide_path_ok or not records[acceptor]["unambiguous"] or
          not records[donor-1]["unambiguous"] or not records[donor]["unambiguous"]):
        hbond_flags.append(False)
      else:
        hbond_flags.append(_odd_backbone_hbond(records, acceptor, donor))
    all_supported = bool(hbond_flags) and all(hbond_flags)
    if all_supported:
      label += " (H-bond supported)"
    supported_features.append((first, last, label, pattern_count, all_supported))
  return supported_features


ODD_HELICAL_FEATURE_COLOUR_INDEX = 28
ODD_HELICAL_HIGHLIGHT_STATE = {}


def _clear_odd_helical_feature_highlights(highlight_state=None):
  highlight_state = ODD_HELICAL_HIGHLIGHT_STATE if highlight_state is None else highlight_state
  for mol_id, state in list(highlight_state.items()):
    try:
      if state.get("handles"):
        for handle in state["handles"]:
          try:
            delete_additional_representation(mol_id, handle)
          except Exception:
            pass
        current_handles = CUSTOM_COLOUR_ADDITIONAL_REPRESENTATIONS.get(mol_id, [])
        CUSTOM_COLOUR_ADDITIONAL_REPRESENTATIONS[mol_id] = [
          handle for handle in current_handles if handle not in state["handles"]]
      elif _supports_direct_user_defined_colouring():
        representation = state.get("representation")
        _restore_last_standard_representation(mol_id, representation)
        blank_rows = state.get("blank_rows")
        if blank_rows is None:
          blank_rows = _direct_user_defined_blank_rows_for_molecule(mol_id)
        if blank_rows:
          _set_user_defined_atom_colour_by_residue_rows_py(mol_id, blank_rows)
        _restore_last_standard_representation(mol_id, representation)
    except Exception:
      traceback.print_exc()
    highlight_state.pop(mol_id, None)


def _highlight_odd_helical_range(navigation, highlight_state=None):
  if globals().get("_USER_DEFINED_COLOUR_JOB_ACTIVE", False):
    add_status_bar_text("Please wait for colouring to finish before highlighting a region.")
    return False
  highlight_state = ODD_HELICAL_HIGHLIGHT_STATE if highlight_state is None else highlight_state
  mol_id = navigation.get("mol_id")
  residue_specs = navigation.get("residue_specs") or []
  if mol_id is None or not residue_specs:
    return False
  existing = highlight_state.get(mol_id)
  if existing and _supports_direct_user_defined_colouring():
    previous_specs = existing.get("residue_specs", [])
    previous_keys = {tuple(spec[:3]) for spec in previous_specs}
    next_keys = {tuple(spec[:3]) for spec in residue_specs}
    changed_rows = [
      (spec, 0) for spec in previous_specs if tuple(spec[:3]) not in next_keys
    ] + [
      (spec, ODD_HELICAL_FEATURE_COLOUR_INDEX)
      for spec in residue_specs if tuple(spec[:3]) not in previous_keys
    ]
    if changed_rows and not _set_user_defined_atom_colour_by_residue_rows_py(mol_id, changed_rows):
      return False
    get_bond_type = getattr(coot, "get_graphics_molecule_bond_type", None)
    try:
      ca_view_active = callable(get_bond_type) and get_bond_type(mol_id) == 13
    except Exception:
      ca_view_active = False
    if changed_rows or not ca_view_active:
      _graphics_to_legacy_user_defined_representation(mol_id, changed_rows, ca_only=True)
    existing["residue_specs"] = list(residue_specs)
    return True
  if existing:
    _clear_odd_helical_feature_highlights(highlight_state)
  rows = [(spec, ODD_HELICAL_FEATURE_COLOUR_INDEX) for spec in residue_specs]
  try:
    representation = _last_standard_representation_function(mol_id)
    if _supports_direct_user_defined_colouring():
      ensure_user_defined_colour_table()
      _clear_custom_colour_additional_representations(mol_id, clear_user_colours=False)
      blank_rows = _direct_user_defined_blank_rows_for_molecule(mol_id)
      result = _set_user_defined_atom_colour_by_residue_rows_py(
        mol_id, blank_rows + rows)
      if result:
        _graphics_to_legacy_user_defined_representation(mol_id, blank_rows + rows, ca_only=True)
    else:
      blank_rows = None
      result = _apply_user_defined_residue_colours(mol_id, [], rows)
    if not result:
      return False
    handles = list(CUSTOM_COLOUR_ADDITIONAL_REPRESENTATIONS.get(mol_id, [])) \
      if not _supports_direct_user_defined_colouring() else []
    highlight_state[mol_id] = {
      "representation": representation,
      "handles": handles,
      "residue_specs": list(residue_specs),
      "blank_rows": blank_rows,
    }
    return True
  except Exception:
    traceback.print_exc()
    return False


def _odd_weak_region_hits(records):
  """Group existing backbone support scores within connected protein stretches."""
  regions = []
  chunks = []
  for index, record in enumerate(records):
    usable = (record.get("backbone_support") is not None and
              record["name"] in _REGISTER_PROTEIN_NAMES and "CA" in record["xyz"])
    if not usable:
      continue
    connected = False
    if chunks and chunks[-1][-1] == index - 1:
      previous = records[index - 1]["xyz"]
      current = record["xyz"]
      connected = ("C" in previous and "N" in current and
                   1.0**2 <= _distance_sq(previous["C"], current["N"]) <= 1.8**2)
    if connected:
      chunks[-1].append(index)
    else:
      chunks.append([index])
  for chunk in chunks:
    weak = [records[index]["backbone_support"] < EMRINGER_HELPER_MIN_BACKBONE_SUPPORT_FRACTION
            for index in chunk]
    spans = []
    for start in range(len(chunk) - 4):
      if sum(weak[start:start + 5]) < 3:
        continue
      end = start + 4
      # Trim stronger window edges instead of inflating region boundaries.
      while not weak[start]:
        start += 1
      while not weak[end]:
        end -= 1
      if spans and start <= spans[-1][1] + 2:
        spans[-1][1] = max(spans[-1][1], end)
      else:
        spans.append([start, end])
    for start, end in spans:
      if end - start + 1 < 5:
        continue
      support = [records[index]["backbone_support"] for index in chunk[start:end + 1]]
      regions.append((chunk[start], chunk[end], sum(weak[start:end + 1]),
                      sum(support) / len(support)))
  return regions


def _chain_residue_count_text(records):
  """Describe a chain from existing snapshot rows, not per-residue API calls."""
  polymer = sum(record["name"].strip().upper() in POLYMER_RESIDUE_NAMES
                for _spec, record in records)
  return f"{polymer} polymer, {len(records) - polymer} non-polymer"


def _iter_odd_residue_dialog_data():
  """Build the current Odd residues dialog payload from the active model/map."""
  map_id = _active_analysis_map_or_status()
  if map_id is None:
    return None

  mol_id = _active_polymer_molecule_or_status()
  if mol_id is None:
    return None

  peak_threshold = _odd_residue_peak_threshold_or_status(map_id)
  if peak_threshold is None:
    return None
  misfit_peak_threshold = peak_threshold * EMRINGER_HELPER_MISFIT_MIN_PEAK_TO_CONTOUR_RATIO
  misfit_gain_threshold = peak_threshold * EMRINGER_HELPER_MISFIT_MIN_GAIN_TO_CONTOUR_RATIO
  occupied_distance_sq = EMRINGER_HELPER_OCCUPIED_DISTANCE * EMRINGER_HELPER_OCCUPIED_DISTANCE

  categorized_entries = {category_name: [] for category_name in ODD_RESIDUE_CATEGORY_ORDER}
  missing_atom_residue_keys = _missing_atom_residue_keys(mol_id)
  snapshot = _model_coordinate_snapshot(mol_id)
  chain_records = {}
  for spec, record in snapshot.items():
    chain_records.setdefault(spec[0], []).append((spec, record))
  edit_epoch = _MODEL_EDIT_EPOCHS.get(mol_id, 0)
  processed = 0

  for chain_index, (chain_id, records) in enumerate(chain_records.items()):
    backbone_records = []
    for serial_number, (spec, record) in enumerate(records):
      if processed % 8 == 0:
        yield (processed, len(snapshot))
        if not valid_model_molecule_qm(mol_id) or not valid_map_molecule_qm(map_id):
          raise RuntimeError("The analysis model or map was closed.")
        if (_MODEL_EDIT_EPOCHS.get(mol_id, 0) != edit_epoch or
            _odd_residue_peak_threshold_or_status(map_id) != peak_threshold):
          raise RuntimeError("Model or contour changed; please restart Odd residues.")
      processed += 1
      residue_name_here = record["name"]
      _chain, resno, ins_code = spec

      atom_records, atom_xyz = _residue_atom_records_and_xyz(mol_id, chain_id, resno, ins_code, snapshot=snapshot)
      backbone_records.append({"name": residue_name_here, "xyz": atom_xyz, "resno": resno,
                               "ins": ins_code, "serial": serial_number,
                               "backbone_support": None,
                               "unambiguous": all(sum(atom["name"] == name for atom in atom_records) == 1
                                                  for name in ("N", "CA", "C", "O"))})
      if not atom_xyz:
        continue
      residue_id = _format_residue_id(chain_id, resno, ins_code)
      residue_point = _residue_display_point(atom_records, atom_xyz)
      dialog_label = _odd_residue_dialog_label(chain_id, resno, ins_code, residue_name_here)
      if (chain_id, resno, ins_code or "") in missing_atom_residue_keys:
        dialog_label = f"{dialog_label} (Missing atoms)"
      navigation_metadata = _odd_residue_navigation_metadata(mol_id, chain_id, resno, ins_code, serial_number)

      def add_result(category_name, sort_key, detail_label, point=None, use_navigation=False):
        _append_odd_residue_entry(
          categorized_entries,
          category_name,
          sort_key,
          dialog_label,
          point if point is not None else residue_point,
          detail_label,
          navigation_metadata if use_navigation else {"chain_id": chain_id},
        )

      if residue_name_here not in POLYMER_RESIDUE_NAMES:
        density_summary = _heavy_atom_density_summary(map_id, atom_records, peak_threshold)
        if density_summary["count"] == 0:
          continue
        if residue_name_here in WATER_RESIDUE_NAMES or density_summary["count"] == 1:
          max_density = density_summary["max_density"]
          if max_density < peak_threshold:
            add_result(
              "Weak Waters/Ions",
              (max_density, chain_index, serial_number),
              f"{residue_id} {residue_name_here}: max atom density {max_density:.4f} "
              f"below threshold {peak_threshold:.4f}",
              point=density_summary["max_atom"]["xyz"],
            )
          continue
        outside_fraction = density_summary["outside_count"] / float(density_summary["count"])
        if outside_fraction > EMRINGER_HELPER_WEAK_LIGAND_OUTSIDE_FRACTION:
          mean_density = density_summary["mean_density"]
          add_result(
            "Weak Ligands",
            (-outside_fraction, mean_density, chain_index, serial_number),
            f"{residue_id} {residue_name_here}: {density_summary['outside_count']}/{density_summary['count']} heavy atoms "
            f"({outside_fraction:.0%}) outside contour; mean density {mean_density:.4f}",
          )
        continue

      if residue_name_here == "GLY" and backbone_records[-1]["unambiguous"]:
        # Keep the L-amino-acid handedness; the carbonyl oxygen should not
        # choose which of glycine's two alpha-hydrogen sites is tested.
        cb_direction = _pseudo_cb_direction_from_backbone(
          {name: atom_xyz[name] for name in ("N", "CA", "C")})
        if cb_direction is not None:
          cb_point = [atom_xyz["CA"][axis] + 1.53 * cb_direction[axis] for axis in range(3)]
          cb_density = density_at_point(map_id, *cb_point)
          if math.isfinite(cb_density) and cb_density >= peak_threshold:
            _append_odd_residue_entry(
              categorized_entries, "Possible Misfits",
              (0.0, -cb_density, chain_index, serial_number, -1),
              f"{dialog_label}: unexpected C-beta density", residue_point,
              f"{residue_id} GLY: unexpected C-beta density; possible extra sidechain density; "
              f"density {cb_density:.4f} at/above threshold {peak_threshold:.4f}. "
              "Check residue identity and nearby atoms; this is a tentative flag.",
              navigation_metadata,
            )

      support_fraction, supported_points, total_points = _backbone_density_support(map_id, atom_xyz, peak_threshold)
      if total_points:
        backbone_records[-1]["backbone_support"] = support_fraction
      if total_points == 0 or support_fraction < EMRINGER_HELPER_MIN_BACKBONE_SUPPORT_FRACTION:
        add_result(
          "Weak Backbone",
          (support_fraction, chain_index, serial_number),
          f"{residue_id} {residue_name_here}: backbone support {supported_points}/{total_points} "
          f"({support_fraction:.0%}) below threshold",
          use_navigation=True,
        )
        continue

      if residue_name_here == "GLY":
        continue

      stage_specs = _emringer_stage_specs(atom_xyz, atom_records, residue_name_here)
      if not stage_specs:
        continue

      stage_scan = _scan_polymer_stages(
        map_id,
        residue_id,
        residue_name_here,
        stage_specs,
        peak_threshold,
        misfit_peak_threshold,
        misfit_gain_threshold,
        occupied_distance_sq,
      )

      for stage_index, stage_result in stage_scan["possible_misfits"]:
        add_result(
          "Possible Misfits",
          (-stage_result["density_gain"], -stage_result["peak"]["density"], chain_index, serial_number, stage_index),
          _possible_misfit_detail_label(residue_id, residue_name_here, stage_result, peak_threshold),
          point=stage_result["peak"]["point"],
          use_navigation=True,
        )

      best_stage_result = stage_scan["best_stage_result"]
      if best_stage_result is None or best_stage_result["peak"]["density"] >= peak_threshold:
        continue

      if not stage_scan["has_stage_atoms"]:
        continue

      strongest_buffered_peak_density = stage_scan["strongest_buffered_peak_density"]
      if strongest_buffered_peak_density is not None and strongest_buffered_peak_density >= peak_threshold:
        continue

      strongest_current_density = stage_scan["strongest_current_density"]
      if strongest_current_density is not None and strongest_current_density >= peak_threshold:
        continue

      if strongest_current_density is None:
        current_density_note = "no current stage atoms placed"
      else:
        current_density_note = "strongest current stage density {0:.4f}".format(strongest_current_density)
      add_result(
        "Weak Sidechains",
        (best_stage_result["peak"]["density"], -(strongest_current_density or 0.0), chain_index, serial_number),
        f"{residue_id} {residue_name_here}: best {best_stage_result['label']} peak "
        f"{best_stage_result['peak']['density']:.4f} below threshold {peak_threshold:.4f}; "
        f"{current_density_note}; strongest inset-ring peak {0:.4f}".format(strongest_buffered_peak_density or 0.0),
        point=(best_stage_result["peak"]["point"] or residue_point),
        use_navigation=True,
      )

    for first, last, weak_count, mean_support in _odd_weak_region_hits(backbone_records):
      start, end = backbone_records[first], backbone_records[last]
      count = last - first + 1
      label = (f"{_format_residue_id(chain_id, start['resno'], start['ins'])} - "
               f"{_format_residue_id(chain_id, end['resno'], end['ins'])} | "
               f"{count} residues | BB support {mean_support:.0%}")
      middle = backbone_records[(first + last) // 2]
      navigation = _odd_residue_navigation_metadata(
        mol_id, chain_id, middle["resno"], middle["ins"], middle["serial"])
      navigation["type"] = "polymer_range"
      navigation["residue_specs"] = [
        [chain_id, record["resno"], record["ins"] or ""]
        for record in backbone_records[first:last + 1]]
      _append_odd_residue_entry(
        categorized_entries, "Weak regions", (mean_support, chain_index, first), label,
        middle["xyz"]["CA"],
        f"{label}; {weak_count}/{count} residues meet the Weak Backbone criterion "
        f"(support below {EMRINGER_HELPER_MIN_BACKBONE_SUPPORT_FRACTION:.0%} at the current contour).",
        navigation)

    for index, context, phi, psi in _odd_backbone_context_hits(backbone_records):
      record = backbone_records[index]
      detail = (f"{_format_residue_id(chain_id, record['resno'], record['ins'])} {record['name']}: "
                f"phi {phi:.1f}, psi {psi:.1f}; possible disruption within a {context}-like backbone run")
      _append_odd_residue_entry(
        categorized_entries, "Unusual phi/psi in context", (chain_index, record["serial"]),
        _odd_residue_dialog_label(chain_id, record["resno"], record["ins"], record["name"]),
        record["xyz"]["CA"], detail,
        _odd_residue_navigation_metadata(mol_id, chain_id, record["resno"], record["ins"], record["serial"]))

    for first, last, kind, pattern_count, hbond_supported in _odd_helical_feature_hits(backbone_records):
      record, end = backbone_records[first], backbone_records[last]
      start_id = _format_residue_id(chain_id, record["resno"], record["ins"])
      end_id = _format_residue_id(chain_id, end["resno"], end["ins"])
      label = f"{start_id} - {end_id}: {kind}"
      confidence_note = "All expected backbone H-bonds also fit the geometry." if hbond_supported else \
                       "No consistent supporting backbone H-bond pattern was found."
      detail = (f"{label}; {pattern_count} consecutive C-alpha geometry windows. "
                f"{confidence_note} Tentative geometry flag, not a validation error.")
      middle_index = (first+last)//2
      middle = backbone_records[middle_index]
      navigation = _odd_residue_navigation_metadata(
        mol_id, chain_id, middle["resno"], middle["ins"], middle["serial"])
      navigation["type"] = "polymer_range"
      navigation["residue_specs"] = [
        [chain_id, item["resno"], item["ins"] or ""] for item in backbone_records[first:last+1]]
      _append_odd_residue_entry(
        categorized_entries, "Helical features", (chain_index, first), label,
        middle["xyz"]["CA"], detail,
        navigation)

  if snapshot != _model_coordinate_snapshot(mol_id):
    raise RuntimeError("Model coordinates changed; please restart Odd residues.")
  sorted_gui_categories, sorted_detail_categories = _sorted_categorized_outputs(
    categorized_entries,
    ODD_RESIDUE_CATEGORY_ORDER,
  )
  total_hits = sum(len(entries) for entries in categorized_entries.values())
  return {
    "mol_id": mol_id,
    "map_id": map_id,
    "peak_threshold": peak_threshold,
    "gui_categories": sorted_gui_categories,
    "chain_choices": [(chain, f"Chain {chain or '(blank)'}: {_chain_residue_count_text(records)}")
                      for chain, records in chain_records.items()],
    "detail_categories": sorted_detail_categories,
    "title": _odd_residue_dialog_title(mol_id, map_id),
    "total_hits": total_hits,
  }


def _collect_odd_residue_dialog_data():
  iterator = _iter_odd_residue_dialog_data()
  while True:
    try:
      next(iterator)
    except StopIteration as result:
      return result.value


def _start_odd_residue_analysis(on_complete):
  """Yield between small residue batches; never run Coot APIs on worker threads."""
  glib = _coot_gui_repository_module("GLib")
  if glib is None or Gtk is None:
    return on_complete(_collect_odd_residue_dialog_data())
  window = Gtk.Window(title="Odd residues")
  window.set_default_size(340, 100)
  box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=8)
  box.set_margin_top(12)
  box.set_margin_bottom(12)
  box.set_margin_start(12)
  box.set_margin_end(12)
  progress = Gtk.Label(label="Analysing model...")
  cancel = Gtk.Button(label="Cancel")
  box.append(progress)
  box.append(cancel)
  window.set_child(box)
  state = {"cancelled": False}
  iterator = _iter_odd_residue_dialog_data()
  def cancel_analysis(*_args):
    state["cancelled"] = True
    return False
  cancel.connect("clicked", cancel_analysis)
  window.connect("close-request", cancel_analysis)
  def step():
    data = None
    if not state["cancelled"]:
      try:
        count, total = next(iterator)
        progress.set_text(f"Analysing residues: {count}/{total}")
        return True
      except StopIteration as result:
        data = result.value
      except Exception as error:
        info_dialog(f"Odd residues analysis stopped: {error}")
    iterator.close()
    window.destroy()
    on_complete(data)
    return False
  window.present()
  glib.idle_add(step)
  return window


def _show_odd_residue_dialog_data(dialog_data):
  if dialog_data is None:
    return None

  if dialog_data["total_hits"] == 0:
    info_dialog(
      "No odd residues found in the current molecule\n"
      f"at the current threshold ({dialog_data['peak_threshold']:.4f})."
    )
    return 0

  def refresh_payload():
    refreshed_data = _collect_odd_residue_dialog_data()
    if refreshed_data is None:
      return None
    _log_odd_residue_results(
      refreshed_data["mol_id"],
      refreshed_data["map_id"],
      refreshed_data["peak_threshold"],
      refreshed_data["detail_categories"],
    )
    return refreshed_data["title"], refreshed_data["gui_categories"], refreshed_data["chain_choices"]

  def refresh_async(done):
    def completed(data):
      if data is not None:
        _log_odd_residue_results(data["mol_id"], data["map_id"], data["peak_threshold"], data["detail_categories"])
      done((data["title"], data["gui_categories"], data["chain_choices"]) if data else None)
    return _start_odd_residue_analysis(completed)
  refresh_payload.start_async = refresh_async

  categorized_interesting_things_gui(
    dialog_data["title"],
    dialog_data["gui_categories"],
    refresh_function=refresh_payload,
    chain_choices=dialog_data["chain_choices"],
  )
  _log_odd_residue_results(
    dialog_data["mol_id"],
    dialog_data["map_id"],
    dialog_data["peak_threshold"],
    dialog_data["detail_categories"],
  )
  return dialog_data["total_hits"]


def find_odd_residues():
  """Find odd residues without blocking the interface for a whole-model scan."""
  return _start_odd_residue_analysis(_show_odd_residue_dialog_data)


def find_emringer_like_sidechain_density_outliers():
  """Backward-compatible entry point for the older menu/helper name."""
  return find_odd_residues()



RESIDUE_ANNOTATION_DEFAULT_AUTHOR = (
  os.environ.get("USER")
  or os.environ.get("USERNAME")
  or ""
)
RESIDUE_ANNOTATION_LAST_AUTHOR = RESIDUE_ANNOTATION_DEFAULT_AUTHOR
RESIDUE_ANNOTATIONS_BY_MOLECULE = {}
RESIDUE_ANNOTATION_CATEGORY_PREFIX = "_cootnote_residue_note."
RESIDUE_ANNOTATION_META_VERSION_TAG = "_cootnote_annotation_meta.version"
RESIDUE_ANNOTATION_META_PROGRAM_TAG = "_cootnote_annotation_meta.program"



def _annotation_now_utc():
  return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _annotation_text_or_empty(value):
  if value in (None, ".", "?"):
    return ""
  return str(value)


def _annotation_optional_text(value):
  text = _annotation_text_or_empty(value).strip()
  if text in ("''", '""'):
    return ""
  return text


def _annotation_safe_int(value, default=None):
  try:
    return int(str(value).strip())
  except Exception:
    return default


def _annotation_label_seq_text(value):
  return _annotation_text_or_empty(value).strip()


def _annotation_encode_text(value):
  text = _annotation_text_or_empty(value)
  encoded = base64.b64encode(text.encode("utf-8"))
  return encoded.decode("ascii")


def _annotation_decode_text(value):
  text = _annotation_text_or_empty(value)
  if not text:
    return ""
  try:
    return base64.b64decode(text.encode("ascii")).decode("utf-8")
  except Exception:
    return text


def _annotation_residue_key(chain_id, resno, ins_code, residue_name_here=""):
  return (str(chain_id), int(resno), str(ins_code or ""), str(residue_name_here or ""))


def _annotation_source_path(mol_id, mmcif_only=False):
  candidates = []
  for candidate in (
    molecule_name_stub_py(mol_id, 0),
    molecule_name_stub_py(mol_id, 1),
    molecule_name(mol_id),
  ):
    if isinstance(candidate, str) and candidate:
      expanded = os.path.abspath(os.path.expanduser(candidate))
      if expanded not in candidates:
        candidates.append(expanded)
  for candidate in candidates:
    if os.path.isfile(candidate):
      if mmcif_only and not candidate.lower().endswith((".cif", ".mmcif")):
        continue
      return candidate
  return None


def _annotation_source_display_name(mol_id):
  source_path = _annotation_source_path(mol_id)
  if source_path:
    return os.path.basename(source_path)
  return None


def _annotation_export_directory():
  desktop_path = os.path.expanduser("~/Desktop")
  if os.path.isdir(desktop_path):
    return desktop_path
  return _coot_trimmings_dir()


def _annotation_entry_count_label(count):
  return f"{count} annotation{'s' if count != 1 else ''}"


def _annotation_entries_for_molecule(mol_id):
  return RESIDUE_ANNOTATIONS_BY_MOLECULE.setdefault(mol_id, [])


def _annotation_default_author():
  return RESIDUE_ANNOTATION_LAST_AUTHOR or RESIDUE_ANNOTATION_DEFAULT_AUTHOR or ""


def _annotation_next_id(mol_id):
  entries = _annotation_entries_for_molecule(mol_id)
  max_id = 0
  for entry in entries:
    entry_id = _annotation_safe_int(entry.get("id"), default=0)
    if entry_id > max_id:
      max_id = entry_id
  return max_id + 1


def _annotation_prepare_note_fields(title_text, author_text, note_text):
  normalized_title = _annotation_optional_text(title_text)
  normalized_author = _annotation_text_or_empty(author_text).strip() or _annotation_default_author()
  normalized_note = _annotation_text_or_empty(note_text).strip()
  if not normalized_note and not normalized_title:
    info_dialog("Please enter a title or note before saving the annotation.")
    return None
  return normalized_title, normalized_author, normalized_note


def _annotation_note_title_placeholder(note_text, max_length=48):
  preview_entry = {"title": "", "note": _annotation_text_or_empty(note_text)}
  placeholder_text = _annotation_entry_display_title(
    preview_entry,
    fallback_from_note=True,
    max_length=max_length,
  )
  if placeholder_text == "Untitled note":
    return ""
  return placeholder_text


def _annotation_target_for_active_residue(expected_mol_id=None):
  residue = active_residue()
  if not residue:
    info_dialog("Select an active residue first.")
    return None
  mol_id, chain_id, resno, ins_code = residue[:4]
  if expected_mol_id is not None and mol_id != expected_mol_id:
    info_dialog(
      f"The active residue is now in molecule #{mol_id}, not the annotation browser molecule #{expected_mol_id}."
    )
    return None
  residue_name_here = residue_name(mol_id, chain_id, resno, ins_code) or ""
  if not residue_name_here:
    info_dialog("Failed to determine the active residue name.")
    return None
  return {
    "mol_id": mol_id,
    "chain_id": str(chain_id),
    "resno": int(resno),
    "ins_code": str(ins_code or ""),
    "residue_name": residue_name_here,
  }


def _new_annotation_entry(target, title_text, author_text, note_text, entry_id):
  entry = {
    "id": int(entry_id),
    "label_comp_id": target["residue_name"],
    "label_asym_id": _annotation_text_or_empty(target.get("label_asym_id")).strip(),
    "label_seq_id": _annotation_label_seq_text(target.get("label_seq_id")),
    "auth_asym_id": target["chain_id"],
    "auth_seq_id": int(target["resno"]),
    "pdbx_PDB_ins_code": str(target["ins_code"] or ""),
    "title": _annotation_optional_text(title_text),
    "author": _annotation_text_or_empty(author_text).strip() or _annotation_default_author(),
    "modified_utc": _annotation_now_utc(),
    "note": _annotation_text_or_empty(note_text).strip(),
  }
  entry["target_kind"] = "region" if target.get("end") else "residue"
  if target.get("end"):
    end = target["end"]
    entry.update({"end_label_comp_id": end["residue_name"],
                  "end_label_asym_id": end.get("label_asym_id", ""),
                  "end_label_seq_id": end.get("label_seq_id", ""),
                  "end_auth_asym_id": end["chain_id"], "end_auth_seq_id": end["resno"],
                  "end_pdbx_PDB_ins_code": end.get("ins_code", "")})
  return entry


def _annotation_normalize_entry(raw_entry, fallback_id=None):
  label_comp_id = _annotation_text_or_empty(raw_entry.get("label_comp_id")).strip() or "?"
  auth_asym_id = _annotation_text_or_empty(raw_entry.get("auth_asym_id")).strip()
  auth_seq_id = _annotation_safe_int(raw_entry.get("auth_seq_id"), default=None)
  if not auth_asym_id or auth_seq_id is None:
    return None
  entry_id = _annotation_safe_int(raw_entry.get("id"), default=fallback_id)
  if entry_id is None:
    return None
  normalized_entry = {
    "id": entry_id,
    "label_comp_id": label_comp_id,
    "label_asym_id": _annotation_text_or_empty(raw_entry.get("label_asym_id")).strip(),
    "label_seq_id": _annotation_label_seq_text(raw_entry.get("label_seq_id")),
    "auth_asym_id": auth_asym_id,
    "auth_seq_id": auth_seq_id,
    "pdbx_PDB_ins_code": _annotation_text_or_empty(raw_entry.get("pdbx_PDB_ins_code")).strip(),
    "title": _annotation_optional_text(raw_entry.get("title")),
    "author": _annotation_text_or_empty(raw_entry.get("author")).strip() or _annotation_default_author(),
    "modified_utc": _annotation_text_or_empty(raw_entry.get("modified_utc")).strip() or _annotation_now_utc(),
    "note": _annotation_text_or_empty(raw_entry.get("note")).strip(),
  }
  if not normalized_entry["note"] and not normalized_entry["title"]:
    return None
  kind = _annotation_text_or_empty(raw_entry.get("target_kind")).strip() or "residue"
  if kind not in ("residue", "region"):
    return None
  normalized_entry["target_kind"] = kind
  if kind == "region":
    for field in RESIDUE_ANNOTATION_RANGE_FIELDS[1:]:
      normalized_entry[field] = _annotation_text_or_empty(raw_entry.get(field)).strip()
    end_number = _annotation_safe_int(raw_entry.get("end_auth_seq_id"))
    if end_number is None or normalized_entry["end_auth_asym_id"] != auth_asym_id:
      return None
    if (normalized_entry["label_asym_id"] and normalized_entry["end_label_asym_id"]
        and normalized_entry["label_asym_id"] != normalized_entry["end_label_asym_id"]):
      return None
    normalized_entry["end_auth_seq_id"] = end_number
  return normalized_entry


def _annotation_entry_group_key(entry):
  key = _annotation_residue_key(entry["auth_asym_id"], entry["auth_seq_id"],
                                entry.get("pdbx_PDB_ins_code", ""), entry.get("label_comp_id", ""))
  if entry.get("target_kind") == "region":
    key += _annotation_residue_key(entry["end_auth_asym_id"], entry["end_auth_seq_id"],
                                   entry.get("end_pdbx_PDB_ins_code", ""), entry.get("end_label_comp_id", ""))
  return key


def _annotation_entry_label(entry):
  if entry.get("target_kind") == "region":
    return (f"{entry['auth_asym_id']}:{entry['auth_seq_id']}{entry.get('pdbx_PDB_ins_code', '')}"
            f"-{entry['end_auth_seq_id']}{entry.get('end_pdbx_PDB_ins_code', '')}")
  return _odd_residue_dialog_label(entry["auth_asym_id"], entry["auth_seq_id"],
                                 entry.get("pdbx_PDB_ins_code", ""), entry.get("label_comp_id", ""))


def _annotation_region_specs(mol_id, entry, all_specs=None):
  if entry.get("target_kind") != "region":
    return [(entry["auth_asym_id"], entry["auth_seq_id"], entry.get("pdbx_PDB_ins_code", ""))]
  specs = [tuple(s) for s in (_all_residue_specs_for_colouring(mol_id) if all_specs is None else all_specs)
           if s[0] == entry["auth_asym_id"]]
  first = (entry["auth_asym_id"], entry["auth_seq_id"], entry.get("pdbx_PDB_ins_code", ""))
  last = (entry["end_auth_asym_id"], entry["end_auth_seq_id"], entry.get("end_pdbx_PDB_ins_code", ""))
  if specs.count(first) != 1 or specs.count(last) != 1:
    return []
  lower, upper = sorted((specs.index(first), specs.index(last)))
  members = [spec for spec in specs[lower:upper+1] if _residue_is_polymer(mol_id, *spec)]
  return members if first in members and last in members else []


def _annotation_entry_display_title(entry, fallback_from_note=True, max_length=None):
  title_text = _annotation_optional_text(entry.get("title")).replace("\r\n", "\n").replace("\r", "\n").strip()
  if title_text:
    title_text = title_text.splitlines()[0].strip()
  elif fallback_from_note:
    note_lines = _annotation_text_or_empty(entry.get("note")).splitlines()
    title_text = note_lines[0].strip() if note_lines else ""
  if max_length is not None and len(title_text) > max_length:
    title_text = title_text[:max_length - 1].rstrip() + "..."
  return title_text


def _annotation_group_sort_key(group):
  first_entry = group["entries"][0]
  return (group["key"][0], group["key"][1], group["key"][2], first_entry.get("id", 0))


def _annotation_entry_sort_key(entry):
  return (entry.get("modified_utc", ""), entry.get("id", 0))


def _annotation_groups_for_molecule(mol_id):
  grouped = {}
  for entry in _annotation_entries_for_molecule(mol_id):
    group_key = _annotation_entry_group_key(entry)
    grouped.setdefault(group_key, []).append(entry)

  groups = []
  for group_key, entries in grouped.items():
    entries.sort(key=_annotation_entry_sort_key)
    groups.append({
      "key": group_key,
      "dialog_label": _annotation_entry_label(entries[0]),
      "entries": entries,
    })
  groups.sort(key=_annotation_group_sort_key)
  return groups


def _annotation_browser_data(mol_id):
  groups = _annotation_groups_for_molecule(mol_id)
  return {
    "mol_id": mol_id,
    "title": _annotation_browser_title(mol_id),
    "groups": groups,
  }


def _annotation_group_detail_text(group):
  lines = [group["dialog_label"], ""]
  for entry_index, entry in enumerate(group["entries"], start=1):
    lines.append(f"Note {entry_index}")
    explicit_title = _annotation_entry_display_title(entry, fallback_from_note=False)
    if explicit_title:
      lines.append(f"Title: {explicit_title}")
    lines.append(f"Author: {_annotation_text_or_empty(entry.get('author')) or 'Unknown'}")
    lines.append(f"Timestamp: {_annotation_text_or_empty(entry.get('modified_utc')) or 'Unknown'}")
    lines.append("")
    lines.append(_annotation_text_or_empty(entry.get("note")) or "")
    if entry_index != len(group["entries"]):
      lines.extend(["", "-" * 40, ""])
  return "\n".join(lines)


def _annotation_group_preview_text(group, max_length=56):
  preview_text = _annotation_entry_display_title(group["entries"][-1], fallback_from_note=True, max_length=max_length)
  note_count = len(group["entries"])
  suffix = "note" if note_count == 1 else "notes"
  if preview_text:
    return f"{group['dialog_label']} ({note_count} {suffix}) - {preview_text}"
  return f"{group['dialog_label']} ({note_count} {suffix})"


def _annotation_group_nearby_preview_text(group, max_length=48):
  preview_text = _annotation_group_preview_text(group, max_length=max_length)
  distance = group.get("distance")
  if distance is None:
    return preview_text
  return f"[{distance:.1f} A] {preview_text}"


def _annotation_clear_box_children(box):
  child = box.get_first_child()
  while child is not None:
    next_child = child.get_next_sibling()
    box.remove(child)
    child = next_child


def _annotation_detail_label(text, selectable=False, margin_bottom=0):
  label = Gtk.Label(label=text)
  label.set_xalign(0.0)
  label.set_wrap(True)
  try:
    label.set_selectable(bool(selectable))
  except Exception:
    pass
  label.set_margin_start(12)
  label.set_margin_end(12)
  if margin_bottom:
    label.set_margin_bottom(margin_bottom)
  return label


def _annotation_set_detail_placeholder(detail_box, placeholder_text):
  _annotation_clear_box_children(detail_box)
  placeholder_label = _annotation_detail_label(placeholder_text, selectable=False)
  placeholder_label.set_margin_top(8)
  detail_box.append(placeholder_label)
  return None


def _annotation_note_expander_label(entry, entry_index, total_entries, max_length=56):
  preview_title = _annotation_entry_display_title(entry, fallback_from_note=True, max_length=max_length)
  if preview_title:
    return f"Note {entry_index} of {total_entries}: {preview_title}"
  return f"Note {entry_index} of {total_entries}"


def _annotation_render_group_detail_widgets(detail_box, group):
  _annotation_clear_box_children(detail_box)
  if not group:
    return _annotation_set_detail_placeholder(detail_box, "Select a residue annotation thread to view its notes.")
  if not hasattr(Gtk, "Expander"):
    fallback_label = _annotation_detail_label(_annotation_group_detail_text(group), selectable=True)
    fallback_label.set_margin_top(8)
    detail_box.append(fallback_label)
    return None

  note_count = len(group["entries"])
  summary_label = _annotation_detail_label(
    f"{group['dialog_label']}\n{note_count} {'note' if note_count == 1 else 'notes'}",
    selectable=False,
    margin_bottom=6,
  )
  summary_label.set_margin_top(6)
  detail_box.append(summary_label)

  separator = Gtk.Separator(orientation=Gtk.Orientation.HORIZONTAL)
  separator.set_margin_start(12)
  separator.set_margin_end(12)
  separator.set_margin_bottom(8)
  detail_box.append(separator)

  for entry_index, entry in enumerate(group["entries"], start=1):
    expander = Gtk.Expander(label=_annotation_note_expander_label(entry, entry_index, note_count))
    expander.set_margin_start(12)
    expander.set_margin_end(12)
    expander.set_margin_bottom(8)
    expander.set_expanded(entry_index == note_count)

    note_box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=4)
    note_box.set_margin_top(6)
    note_box.set_margin_bottom(6)

    explicit_title = _annotation_entry_display_title(entry, fallback_from_note=False)
    if explicit_title:
      note_box.append(_annotation_detail_label(f"Title: {explicit_title}", selectable=True))
    note_box.append(
      _annotation_detail_label(
        f"Author: {_annotation_text_or_empty(entry.get('author')) or 'Unknown'}",
        selectable=True,
      )
    )
    note_box.append(
      _annotation_detail_label(
        f"Timestamp: {_annotation_text_or_empty(entry.get('modified_utc')) or 'Unknown'}",
        selectable=True,
        margin_bottom=4,
      )
    )
    if _annotation_text_or_empty(entry.get("note")).strip():
      note_box.append(_annotation_detail_label("Note:", selectable=False))
      note_text_label = _annotation_detail_label(_annotation_text_or_empty(entry.get("note")), selectable=True)
      note_text_label.set_margin_start(24)
      note_box.append(note_text_label)

    expander.set_child(note_box)
    detail_box.append(expander)
  return None


def _annotation_export_entry_maps(entity_records):
  polymer_residue_map = {}
  nonpoly_residue_map = {}
  for record in entity_records or []:
    if record["kind"] == "polymer":
      polymer_residue_map[record["chain_id"]] = record
    elif record["kind"] == "nonpolymer":
      nonpoly_residue_map[record["residue_key"]] = record
  return polymer_residue_map, nonpoly_residue_map


def _annotation_export_label_fields(entry, polymer_residue_map=None, nonpoly_residue_map=None):
  label_comp_id = _annotation_text_or_empty(entry.get("label_comp_id")).strip() or "?"
  label_asym_id = _annotation_text_or_empty(entry.get("label_asym_id")).strip()
  label_seq_id = _annotation_label_seq_text(entry.get("label_seq_id"))
  if polymer_residue_map is None or nonpoly_residue_map is None:
    return label_comp_id, label_asym_id, label_seq_id
  residue_key = (
    _annotation_text_or_empty(entry.get("auth_asym_id")).strip(),
    str(entry.get("auth_seq_id", "")),
    _annotation_text_or_empty(entry.get("pdbx_PDB_ins_code")).strip(),
    label_comp_id,
  )
  nonpoly_record = nonpoly_residue_map.get(residue_key)
  if nonpoly_record is not None:
    return label_comp_id, nonpoly_record["asym_id"], "."

  polymer_record = polymer_residue_map.get(residue_key[0])
  if polymer_record is not None:
    exported_label_seq_id = polymer_record["residue_key_to_label_seq"].get(residue_key)
    if exported_label_seq_id is not None:
      return label_comp_id, polymer_record["asym_id"], str(exported_label_seq_id)

  return label_comp_id, label_asym_id, label_seq_id


def _annotation_mmcif_rows(mol_id, entity_records=None):
  rows = []
  polymer_residue_map, nonpoly_residue_map = _annotation_export_entry_maps(entity_records)
  for entry in sorted(_annotation_entries_for_molecule(mol_id), key=_annotation_entry_sort_key):
    label_comp_id, label_asym_id, label_seq_id = _annotation_export_label_fields(
      entry,
      polymer_residue_map=polymer_residue_map,
      nonpoly_residue_map=nonpoly_residue_map,
    )
    rows.append(
      [
        _annotation_text_or_empty(entry.get("id")) or "?",
        label_comp_id or "?",
        label_asym_id or (_annotation_text_or_empty(entry.get("auth_asym_id")) or "?"),
        label_seq_id or ".",
        _annotation_text_or_empty(entry.get("auth_asym_id")) or "?",
        str(entry.get("auth_seq_id", "")),
        _annotation_text_or_empty(entry.get("pdbx_PDB_ins_code")) or "?",
        _annotation_text_or_empty(entry.get("title")),
        _annotation_text_or_empty(entry.get("author")),
        _annotation_text_or_empty(entry.get("modified_utc")) or "?",
        _annotation_text_or_empty(entry.get("note")),
      ]
    )
    end_entry = {field: entry.get("end_" + field) for field in
                 ("label_comp_id", "label_asym_id", "label_seq_id", "auth_asym_id",
                  "auth_seq_id", "pdbx_PDB_ins_code")}
    end_comp, end_asym, end_seq = _annotation_export_label_fields(
      end_entry, polymer_residue_map, nonpoly_residue_map)
    rows[-1].extend([entry.get("target_kind", "residue"),
                    end_comp if entry.get("target_kind") == "region" else "?",
                    end_asym or "?", end_seq or ".",
                    _annotation_text_or_empty(entry.get("end_auth_asym_id")) or "?",
                    _annotation_text_or_empty(entry.get("end_auth_seq_id")) or "?",
                    _annotation_text_or_empty(entry.get("end_pdbx_PDB_ins_code")) or "?"])
  return rows


def _annotation_markdown_cell(value):
  text = _annotation_text_or_empty(value).strip()
  if not text:
    return ""
  return (
    text
    .replace("\\", "\\\\")
    .replace("|", "\\|")
    .replace("\r\n", "\n")
    .replace("\r", "\n")
    .replace("\n", "<br>")
  )


def _annotation_markdown_table_text(mol_id):
  source_display_name = _annotation_source_display_name(mol_id) or "Unknown source"
  lines = [
    f"# Residue annotations for molecule #{mol_id}",
    "",
    f"Source: `{source_display_name}`",
    "",
    "| Residue | Title | Author | Timestamp | Note |",
    "| --- | --- | --- | --- | --- |",
  ]
  for group in _annotation_groups_for_molecule(mol_id):
    residue_label = _annotation_markdown_cell(group["dialog_label"])
    for entry in group["entries"]:
      lines.append(
        "| {0} | {1} | {2} | {3} | {4} |".format(
          residue_label,
          _annotation_markdown_cell(entry.get("title")) or "",
          _annotation_markdown_cell(entry.get("author")) or "Unknown",
          _annotation_markdown_cell(entry.get("modified_utc")) or "",
          _annotation_markdown_cell(entry.get("note")) or "",
        )
      )
  return "\n".join(lines) + "\n"


def _annotation_clear_existing_mmcif_items(block):
  loop_item = block.find_loop_item(RESIDUE_ANNOTATION_CATEGORY_PREFIX + "id")
  if loop_item is not None:
    loop_item.erase()
  for tag in (RESIDUE_ANNOTATION_META_VERSION_TAG, RESIDUE_ANNOTATION_META_PROGRAM_TAG):
    pair_item = block.find_pair_item(tag)
    if pair_item is not None:
      pair_item.erase()


def _load_annotations_from_mmcif_file(file_path):
  if gemmi is None:
    raise RuntimeError("gemmi is not available in this Coot Python environment.")
  document = gemmi.cif.read_file(file_path)
  block = document.sole_block()
  column = block.find_loop(RESIDUE_ANNOTATION_CATEGORY_PREFIX + "id")
  if column is None:
    return []
  loop = column.get_loop()
  if loop is None:
    return []

  tag_indices = {tag: index for index, tag in enumerate(loop.tags)}
  required_tags = [
    RESIDUE_ANNOTATION_CATEGORY_PREFIX + "id",
    RESIDUE_ANNOTATION_CATEGORY_PREFIX + "label_comp_id",
    RESIDUE_ANNOTATION_CATEGORY_PREFIX + "auth_asym_id",
    RESIDUE_ANNOTATION_CATEGORY_PREFIX + "auth_seq_id",
    RESIDUE_ANNOTATION_CATEGORY_PREFIX + "pdbx_PDB_ins_code",
    RESIDUE_ANNOTATION_CATEGORY_PREFIX + "modified_utc",
  ]
  for required_tag in required_tags:
    if required_tag not in tag_indices:
      raise RuntimeError(f"Embedded annotation loop is missing required tag {required_tag}.")

  def loop_text_field(row_index, plain_tag_name, legacy_b64_tag_name=None):
    if plain_tag_name in tag_indices:
      return _annotation_text_or_empty(gemmi.cif.as_string(loop[row_index, tag_indices[plain_tag_name]]))
    if legacy_b64_tag_name is not None and legacy_b64_tag_name in tag_indices:
      return _annotation_decode_text(loop[row_index, tag_indices[legacy_b64_tag_name]])
    return ""

  author_tag = None
  for candidate_tag in (
    RESIDUE_ANNOTATION_CATEGORY_PREFIX + "author",
    RESIDUE_ANNOTATION_CATEGORY_PREFIX + "author_b64",
  ):
    if candidate_tag in tag_indices:
      author_tag = candidate_tag
      break
  title_tag = None
  for candidate_tag in (
    RESIDUE_ANNOTATION_CATEGORY_PREFIX + "title",
    RESIDUE_ANNOTATION_CATEGORY_PREFIX + "title_b64",
  ):
    if candidate_tag in tag_indices:
      title_tag = candidate_tag
      break
  note_tag = None
  for candidate_tag in (
    RESIDUE_ANNOTATION_CATEGORY_PREFIX + "note",
    RESIDUE_ANNOTATION_CATEGORY_PREFIX + "note_b64",
  ):
    if candidate_tag in tag_indices:
      note_tag = candidate_tag
      break
  if author_tag is None or note_tag is None:
    raise RuntimeError("Embedded annotation loop is missing author or note text fields.")

  entries = []
  for row_index in range(loop.length()):
    identity = {name: loop_text_field(row_index, RESIDUE_ANNOTATION_CATEGORY_PREFIX + name)
                for name in ("id", "label_comp_id", "label_asym_id", "label_seq_id",
                             "auth_asym_id", "auth_seq_id", "pdbx_PDB_ins_code", "modified_utc",
                             *RESIDUE_ANNOTATION_RANGE_FIELDS)}
    normalized_entry = _annotation_normalize_entry(
      {
        **identity,
        "title": loop_text_field(
          row_index,
          RESIDUE_ANNOTATION_CATEGORY_PREFIX + "title",
          RESIDUE_ANNOTATION_CATEGORY_PREFIX + "title_b64",
        ) if title_tag is not None else "",
        "author": loop_text_field(
          row_index,
          RESIDUE_ANNOTATION_CATEGORY_PREFIX + "author",
          RESIDUE_ANNOTATION_CATEGORY_PREFIX + "author_b64",
        ),
        "note": loop_text_field(
          row_index,
          RESIDUE_ANNOTATION_CATEGORY_PREFIX + "note",
          RESIDUE_ANNOTATION_CATEGORY_PREFIX + "note_b64",
        ),
      },
      fallback_id=row_index + 1,
    )
    if normalized_entry is None:
      print(f"Skipping malformed residue annotation row {row_index + 1} in {file_path}")
      continue
    entries.append(normalized_entry)
  return entries


def _annotation_write_entries_to_mmcif_block(block, entries):
  _annotation_clear_existing_mmcif_items(block)
  block.set_pair(RESIDUE_ANNOTATION_META_VERSION_TAG, RESIDUE_ANNOTATION_SCHEMA_VERSION)
  block.set_pair(RESIDUE_ANNOTATION_META_PROGRAM_TAG, "coot-trimmings")
  if entries:
    loop = block.init_loop(
      RESIDUE_ANNOTATION_CATEGORY_PREFIX,
      [
        "id",
        "label_comp_id",
        "label_asym_id",
        "label_seq_id",
        "auth_asym_id",
        "auth_seq_id",
        "pdbx_PDB_ins_code",
        "title",
        "author",
        "modified_utc",
        "note",
      ] + list(RESIDUE_ANNOTATION_RANGE_FIELDS),
    )
    for row in entries:
      if len(row) == 11:
        row = list(row) + ["residue", "?", "?", ".", "?", "?", "?"]
      # Gemmi loop values are CIF tokens, not unquoted application strings.
      loop.add_row([value if index in (3, 6, 12, 13, 14, 15, 16, 17) and value in (".", "?")
                    else gemmi.cif.quote(str(value)) for index, value in enumerate(row)])


def _write_annotations_into_mmcif_file(file_path, entries):
  if gemmi is None:
    raise RuntimeError("gemmi is not available in this Coot Python environment.")
  document = gemmi.cif.read_file(file_path)
  _annotation_write_entries_to_mmcif_block(document.sole_block(), entries)
  document.write_file(file_path)


def _annotation_atom_invalid_for_export(atom):
  atom_name = (getattr(atom, "name", "") or "").strip()
  element_name = getattr(getattr(atom, "element", None), "name", "") or ""
  return (
    not atom_name
    or atom_name in (".", "?")
    or element_name in ("", ".", "?", "X")
  )


def _annotation_clear_loop_category(block, first_tag):
  loop_item = block.find_loop_item(first_tag)
  if loop_item is not None:
    loop_item.erase()


def _annotation_prepare_mmcif_table(block, prefix, tags):
  expected_tags = [prefix + tag for tag in tags]
  table = block.find_mmcif_category(prefix)
  if table is not None and table.width() == len(tags) and list(table.tags) == expected_tags:
    while len(table):
      table.remove_row(len(table) - 1)
    return table
  _annotation_clear_loop_category(block, expected_tags[0])
  pair_item = block.find_pair_item(expected_tags[0])
  if pair_item is not None:
    pair_item.erase()
  block.init_loop(prefix, tags)
  return block.find_mmcif_category(prefix)


def _annotation_polymer_type_cif_name(polymer_type):
  polymer_name = getattr(polymer_type, "name", "")
  return {
    "PeptideL": "polypeptide(L)",
    "PeptideD": "polypeptide(D)",
    "Dna": "polydeoxyribonucleotide",
    "Rna": "polyribonucleotide",
    "DnaRnaHybrid": "polydeoxyribonucleotide/polyribonucleotide hybrid",
    "Pna": "peptide nucleic acid",
  }.get(polymer_name, "other")


def _annotation_export_ins_code(seqid):
  ins_code = getattr(seqid, "icode", "") or ""
  if ins_code in (" ", "\x00"):
    return ""
  return ins_code


def _annotation_nonpoly_asym_id(chain_id, record_index):
  chain_prefix = chain_id or "X"
  return f"{chain_prefix}_np_{record_index}"


def _annotation_residue_export_key(chain_id, residue):
  return (
    chain_id or "",
    str(residue.seqid.num),
    _annotation_export_ins_code(residue.seqid),
    residue.name or "",
  )


def _annotation_export_entity_records(structure):
  if len(structure) == 0:
    return []

  records = []
  next_entity_id = 1
  for chain in structure[0]:
    chain_id = chain.name or ""
    if not chain_id:
      continue

    polymer_span = chain.get_polymer()
    polymer_keys = set()
    if len(polymer_span):
      polymer_keys = {_annotation_residue_export_key(chain_id, residue) for residue in polymer_span}
      residue_names = [residue.name for residue in polymer_span]
      records.append({
        "kind": "polymer",
        "chain_id": chain_id,
        "asym_id": chain_id,
        "entity_id": str(next_entity_id),
        "entity_type": "polymer",
        "polymer_type": _annotation_polymer_type_cif_name(polymer_span.check_polymer_type()),
        "one_letter_sequence": polymer_span.make_one_letter_sequence() or "?",
        "sequence": residue_names,
        "residue_key_to_label_seq": {
          _annotation_residue_export_key(chain_id, residue): str(index)
          for index, residue in enumerate(polymer_span, start=1)
        },
      })
      next_entity_id += 1

    nonpoly_record_index = 1
    for residue in chain:
      residue_key = _annotation_residue_export_key(chain_id, residue)
      if residue_key in polymer_keys:
        continue
      records.append({
        "kind": "nonpolymer",
        "chain_id": chain_id,
        "asym_id": _annotation_nonpoly_asym_id(chain_id, nonpoly_record_index),
        "entity_id": str(next_entity_id),
        "entity_type": "water" if residue.name in WATER_RESIDUE_NAMES else "non-polymer",
        "polymer_type": None,
        "one_letter_sequence": None,
        "sequence": [],
        "comp_id": residue.name or "?",
        "residue_key": residue_key,
      })
      next_entity_id += 1
      nonpoly_record_index += 1
  return records


def _annotation_apply_export_entity_metadata(block, structure):
  entity_records = _annotation_export_entity_records(structure)
  if not entity_records:
    return []

  entity_table = _annotation_prepare_mmcif_table(block, "_entity.", ["id", "type"])
  struct_asym_table = _annotation_prepare_mmcif_table(block, "_struct_asym.", ["id", "entity_id"])
  polymer_records = [record for record in entity_records if record["kind"] == "polymer"]
  nonpoly_records = [record for record in entity_records if record["kind"] == "nonpolymer"]
  entity_poly_table = None
  entity_poly_seq_table = None
  if polymer_records:
    entity_poly_table = _annotation_prepare_mmcif_table(
      block,
      "_entity_poly.",
      ["entity_id", "type", "pdbx_strand_id", "pdbx_seq_one_letter_code"],
    )
    entity_poly_seq_table = _annotation_prepare_mmcif_table(block, "_entity_poly_seq.", ["entity_id", "num", "mon_id"])
  nonpoly_table = None
  if nonpoly_records:
    nonpoly_table = _annotation_prepare_mmcif_table(block, "_pdbx_entity_nonpoly.", ["entity_id", "comp_id"])

  for record in entity_records:
    quote = gemmi.cif.quote
    entity_table.append_row([quote(record["entity_id"]), quote(record["entity_type"])])
    struct_asym_table.append_row([quote(record["asym_id"]), quote(record["entity_id"])])
    if record["kind"] == "polymer" and entity_poly_table is not None and entity_poly_seq_table is not None:
      entity_poly_table.append_row([quote(value) if value != "?" else "?" for value in [
        record["entity_id"],
        record["polymer_type"] or "other",
        record["asym_id"],
        record["one_letter_sequence"] or "?",
      ]])
      for index, mon_id in enumerate(record["sequence"], start=1):
        entity_poly_seq_table.append_row([quote(record["entity_id"]), str(index), quote(mon_id) if mon_id else "?"])
    elif record["kind"] == "nonpolymer" and nonpoly_table is not None:
      nonpoly_table.append_row([quote(record["entity_id"]), quote(record["comp_id"]) if record["comp_id"] else "?"])

  return entity_records


def _annotation_apply_export_atom_site_metadata(block, entity_records):
  required_columns = {
    tag: block.find_values(tag)
    for tag in (
      "_atom_site.group_PDB",
      "_atom_site.label_comp_id",
      "_atom_site.label_asym_id",
      "_atom_site.label_entity_id",
      "_atom_site.auth_asym_id",
      "_atom_site.auth_seq_id",
      "_atom_site.pdbx_PDB_ins_code",
      "_atom_site.label_seq_id",
    )
  }
  if any(column is None for column in required_columns.values()):
    return

  group_column = required_columns["_atom_site.group_PDB"]
  comp_id_column = required_columns["_atom_site.label_comp_id"]
  chain_column = required_columns["_atom_site.label_asym_id"]
  entity_column = required_columns["_atom_site.label_entity_id"]
  auth_chain_column = required_columns["_atom_site.auth_asym_id"]
  auth_seq_column = required_columns["_atom_site.auth_seq_id"]
  ins_code_column = required_columns["_atom_site.pdbx_PDB_ins_code"]
  label_seq_column = required_columns["_atom_site.label_seq_id"]
  if not (
    len(group_column)
    == len(comp_id_column)
    == len(chain_column)
    == len(entity_column)
    == len(auth_chain_column)
    == len(auth_seq_column)
    == len(ins_code_column)
    == len(label_seq_column)
  ):
    return

  polymer_residue_map = {}
  nonpoly_residue_map = {}
  for record in entity_records:
    if record["kind"] == "polymer":
      polymer_residue_map[record["chain_id"]] = record
    elif record["kind"] == "nonpolymer":
      nonpoly_residue_map[record["residue_key"]] = record

  label_remap = {}
  auth_remap = {}
  for row_index in range(len(group_column)):
    auth_chain_id = gemmi.cif.as_string(auth_chain_column[row_index])
    ins_code = gemmi.cif.as_string(ins_code_column[row_index])
    residue_key = (
      auth_chain_id,
      gemmi.cif.as_string(auth_seq_column[row_index]),
      "" if ins_code in (".", "?") else ins_code,
      gemmi.cif.as_string(comp_id_column[row_index]),
    )
    old_key = (gemmi.cif.as_string(chain_column[row_index]),
               _annotation_label_seq_text(gemmi.cif.as_string(label_seq_column[row_index])), residue_key[3], residue_key[2])
    nonpoly_record = nonpoly_residue_map.get(residue_key)
    if nonpoly_record is not None:
      chain_column[row_index] = gemmi.cif.quote(nonpoly_record["asym_id"])
      entity_column[row_index] = nonpoly_record["entity_id"]
      label_seq_column[row_index] = "."
    else:
      polymer_record = polymer_residue_map.get(auth_chain_id)
      if polymer_record is None:
        label_seq_column[row_index] = "."
      else:
        chain_column[row_index] = gemmi.cif.quote(polymer_record["asym_id"])
        entity_column[row_index] = polymer_record["entity_id"]
        label_seq_column[row_index] = polymer_record["residue_key_to_label_seq"].get(residue_key, ".")
    new_labels = (chain_column[row_index], label_seq_column[row_index])
    label_remap.setdefault(old_key, set()).add(new_labels)
    auth_remap.setdefault(residue_key, set()).add(new_labels)

  _annotation_remap_connection_labels(block, label_remap, auth_remap)


def _annotation_remap_connection_labels(block, label_remap, auth_remap):
  """Keep link endpoints synchronized with rewritten atom-site identifiers."""
  for partner in ("ptnr1", "ptnr2"):
    prefix = "_struct_conn." + partner + "_"
    asym = block.find_values(prefix + "label_asym_id")
    seq = block.find_values(prefix + "label_seq_id")
    comp = block.find_values(prefix + "label_comp_id")
    auth_asym = block.find_values(prefix + "auth_asym_id")
    auth_seq = block.find_values(prefix + "auth_seq_id")
    ins = block.find_values("_struct_conn.pdbx_" + partner + "_PDB_ins_code")
    for index in range(len(asym)):
      ins_code = gemmi.cif.as_string(ins[index]) if len(ins) else ""
      name = gemmi.cif.as_string(comp[index]) if len(comp) else "?"
      matches = set()
      if len(auth_asym) and len(auth_seq):
        matches = auth_remap.get((gemmi.cif.as_string(auth_asym[index]), gemmi.cif.as_string(auth_seq[index]), ins_code, name), set())
      if not matches:
        matches = label_remap.get((gemmi.cif.as_string(asym[index]),
                  _annotation_label_seq_text(gemmi.cif.as_string(seq[index])) if len(seq) else ".", name, ins_code), set())
      if len(matches) != 1:
        raise ValueError("Cannot unambiguously remap an exported link endpoint.")
      new_asym, new_seq = next(iter(matches))
      asym[index] = new_asym
      if len(seq):
        seq[index] = new_seq


def _annotation_sanitize_export_structure(structure):
  for model in structure:
    for chain in model:
      for residue in chain:
        invalid_atoms = [
          (atom.name, atom.altloc, atom.element)
          for atom in residue
          if _annotation_atom_invalid_for_export(atom)
        ]
        for atom_name, altloc, element in invalid_atoms:
          residue.remove_atom(atom_name, altloc, element)
  try:
    structure.setup_entities()
  except Exception:
    pass
  try:
    structure.add_entity_types(True)
  except Exception:
    pass
  try:
    structure.add_entity_ids(True)
  except Exception:
    pass
  try:
    structure.assign_label_seq_id(True)
  except Exception:
    pass
  return structure


def _annotation_export_document_from_coot_cif(file_path):
  if gemmi is None:
    raise RuntimeError("gemmi is not available in this Coot Python environment.")
  structure = gemmi.read_structure(file_path)
  _annotation_sanitize_export_structure(structure)
  document = structure.make_mmcif_document()
  block = document.sole_block()
  entity_records = _annotation_apply_export_entity_metadata(block, structure)
  _annotation_apply_export_atom_site_metadata(block, entity_records)
  return document, entity_records


def _load_residue_annotations_for_molecule(mol_id, file_path=None, show_status=True):
  annotation_path = file_path or _annotation_source_path(mol_id, mmcif_only=True)
  if not annotation_path:
    if show_status:
      info_dialog("The active molecule does not have an accessible mmCIF source file.")
    return None
  try:
    loaded_entries = _load_annotations_from_mmcif_file(annotation_path)
  except Exception as error:
    if show_status:
      info_dialog(f"Failed to load embedded residue annotations.\n{error}")
    traceback.print_exc()
    return None
  RESIDUE_ANNOTATIONS_BY_MOLECULE[mol_id] = loaded_entries
  if show_status:
    add_status_bar_text(
      f"Loaded {_annotation_entry_count_label(len(loaded_entries))} from {os.path.basename(annotation_path)}"
    )
  return loaded_entries


def _write_model_with_embedded_annotations(mol_id, file_path):
  output_directory = os.path.dirname(os.path.abspath(file_path))
  os.makedirs(output_directory, exist_ok=True)
  source_handle = None
  source_path = None
  export_handle = None
  export_path = None
  try:
    source_handle, source_path = tempfile.mkstemp(suffix=".cif", dir=output_directory)
    os.close(source_handle)
    source_handle = None
    write_cif_file(mol_id, source_path)
    if not os.path.isfile(source_path) or os.path.getsize(source_path) == 0:
      raise RuntimeError("Coot failed to write a temporary mmCIF file.")
    export_document, entity_records = _annotation_export_document_from_coot_cif(source_path)
    _annotation_write_entries_to_mmcif_block(export_document.sole_block(), _annotation_mmcif_rows(mol_id, entity_records))
    export_handle, export_path = tempfile.mkstemp(suffix=".cif", dir=output_directory)
    os.close(export_handle)
    export_handle = None
    export_document.write_file(export_path)
    gemmi.cif.read_file(export_path).sole_block()
    _load_annotations_from_mmcif_file(export_path)
    os.replace(export_path, file_path)
    return file_path
  finally:
    if source_handle is not None:
      os.close(source_handle)
    if export_handle is not None:
      os.close(export_handle)
    if source_path and os.path.exists(source_path):
      try:
        os.remove(source_path)
      except Exception:
        pass
    if export_path and os.path.exists(export_path):
      try:
        os.remove(export_path)
      except Exception:
        pass


def _default_annotation_export_path(mol_id):
  export_directory = _annotation_export_directory()
  source_path = _annotation_source_path(mol_id)
  if source_path:
    source_stem = os.path.splitext(os.path.basename(source_path))[0]
    return os.path.join(export_directory, f"{source_stem}.annotated.cif")
  source_stub = molecule_name_stub_py(mol_id, 0) or f"molecule_{mol_id}"
  return os.path.join(export_directory, f"{source_stub}.annotated.cif")


def _default_annotation_markdown_export_path(mol_id):
  export_directory = _annotation_export_directory()
  source_path = _annotation_source_path(mol_id)
  if source_path:
    source_stem = os.path.splitext(os.path.basename(source_path))[0]
    return os.path.join(export_directory, f"{source_stem}.annotations.md")
  source_stub = molecule_name_stub_py(mol_id, 0) or f"molecule_{mol_id}"
  return os.path.join(export_directory, f"{source_stub}.annotations.md")


def _annotation_normalize_export_path(raw_output_path, allowed_suffixes, default_suffix, path_label):
  raw_path = str(raw_output_path or "").strip()
  if not raw_path:
    info_dialog(f"Please provide an output {path_label} path.")
    return None
  normalized_path = os.path.abspath(os.path.expanduser(raw_path))
  if os.path.isdir(normalized_path):
    info_dialog(f"Please choose a file name, not a directory, for the {path_label} export.")
    return None
  if not normalized_path.lower().endswith(tuple(allowed_suffixes)):
    normalized_path += default_suffix
  return normalized_path


def _annotation_confirm_export_overwrite(dialog_title, normalized_path, retry_function):
  generic_confirm_dialog(
    dialog_title,
    f"{os.path.basename(normalized_path)} already exists.\nOverwrite it?",
    "Cancel",
    lambda: None,
    "Overwrite",
    retry_function,
  )
  return 0


def _export_residue_annotations_for_molecule(mol_id, output_path, allow_overwrite=False):
  normalized_path = _annotation_normalize_export_path(
    output_path,
    (".cif", ".mmcif"),
    ".cif",
    "mmCIF",
  )
  if normalized_path is None:
    return 0
  if os.path.lexists(normalized_path) and not allow_overwrite:
    return _annotation_confirm_export_overwrite(
      "Overwrite annotated mmCIF",
      normalized_path,
      lambda: _export_residue_annotations_for_molecule(mol_id, normalized_path, allow_overwrite=True),
    )
  try:
    _write_model_with_embedded_annotations(mol_id, normalized_path)
  except Exception as error:
    info_dialog(f"Failed to export annotated mmCIF.\n{error}")
    traceback.print_exc()
    return 0
  add_status_bar_text(
    f"Exported annotated mmCIF with {_annotation_entry_count_label(len(_annotation_entries_for_molecule(mol_id)))} "
    f"to {os.path.basename(normalized_path)}"
  )
  return 1


def _export_residue_annotations_markdown_for_molecule(mol_id, output_path, allow_overwrite=False):
  normalized_path = _annotation_normalize_export_path(
    output_path,
    (".md",),
    ".md",
    "Markdown",
  )
  if normalized_path is None:
    return 0
  if os.path.lexists(normalized_path) and not allow_overwrite:
    return _annotation_confirm_export_overwrite(
      "Overwrite annotation Markdown",
      normalized_path,
      lambda: _export_residue_annotations_markdown_for_molecule(mol_id, normalized_path, allow_overwrite=True),
    )
  try:
    output_directory = os.path.dirname(normalized_path)
    if output_directory:
      os.makedirs(output_directory, exist_ok=True)
    with open(normalized_path, "w", encoding="utf-8") as handle:
      handle.write(_annotation_markdown_table_text(mol_id))
  except Exception as error:
    info_dialog(f"Failed to export annotation Markdown.\n{error}")
    traceback.print_exc()
    return 0
  add_status_bar_text(
    f"Exported annotation Markdown with {_annotation_entry_count_label(len(_annotation_entries_for_molecule(mol_id)))} "
    f"to {os.path.basename(normalized_path)}"
  )
  return 1


def _export_residue_annotations_with_confirmation(mol_id, output_path):
  return _export_residue_annotations_for_molecule(mol_id, output_path, allow_overwrite=False)


def _export_residue_annotations_markdown_with_confirmation(mol_id, output_path):
  return _export_residue_annotations_markdown_for_molecule(mol_id, output_path, allow_overwrite=False)


def _ensure_residue_annotations_loaded(mol_id):
  if mol_id in RESIDUE_ANNOTATIONS_BY_MOLECULE:
    return _annotation_entries_for_molecule(mol_id)
  return _load_residue_annotations_for_molecule(mol_id, show_status=False) or _annotation_entries_for_molecule(mol_id)


def _annotation_active_molecule():
  model_molecules = model_molecule_list()
  residue = active_residue()
  if residue and residue[0] in model_molecules:
    return residue[0]
  mol_id = go_to_atom_molecule_number()
  if mol_id in model_molecules:
    return mol_id
  if model_molecules:
    return model_molecules[0]
  return None


def _annotation_active_molecule_or_status():
  mol_id = _annotation_active_molecule()
  if mol_id is not None:
    return mol_id
  info_dialog("You need an active model molecule for residue annotations.")
  return None


def _append_annotation_for_target(target, title_text, author_text, note_text, refresh_function=None):
  global RESIDUE_ANNOTATION_LAST_AUTHOR
  if not target:
    return 0
  note_fields = _annotation_prepare_note_fields(title_text, author_text, note_text)
  if note_fields is None:
    return 0
  normalized_title, normalized_author, normalized_note = note_fields
  RESIDUE_ANNOTATION_LAST_AUTHOR = normalized_author

  mol_id = target["mol_id"]
  if mol_id not in model_molecule_list():
    info_dialog("The annotation model is no longer available.")
    return 0
  annotation_entry = _new_annotation_entry(
    target,
    normalized_title,
    normalized_author,
    normalized_note,
    _annotation_next_id(mol_id),
  )
  if annotation_entry["target_kind"] == "region":
    if not _annotation_region_specs(mol_id, annotation_entry):
      info_dialog("The region endpoints are no longer available or are ambiguous.")
      return 0
  _annotation_entries_for_molecule(mol_id).append(annotation_entry)
  add_status_bar_text(f"Added annotation for {_annotation_entry_label(annotation_entry)}")
  if callable(refresh_function):
    refresh_function(_annotation_entry_group_key(annotation_entry))
  return 1


def prompt_add_annotation_for_active_residue(mol_id=None, refresh_function=None):
  target = _annotation_target_for_active_residue(expected_mol_id=mol_id)
  if target is None:
    return 0
  residue_label = _odd_residue_dialog_label(
    target["chain_id"],
    target["resno"],
    target["ins_code"],
    target["residue_name"],
  )
  generic_title_author_note_entry(
    "Add residue annotation",
    residue_label,
    "",
    _annotation_default_author(),
    "",
    "Add note",
    lambda title_text, author_text, note_text: _append_annotation_for_target(
      target,
      title_text,
      author_text,
      note_text,
      refresh_function,
    ),
  )
  return 1


def _annotation_parse_range(text):
  match = re.fullmatch(r"\s*(-?\d+)([A-Za-z]*)\s*(?:-|,)\s*(-?\d+)([A-Za-z]*)\s*", str(text))
  if not match:
    raise ValueError("Use a range such as 75-80 or 75A-80B.")
  return (int(match[1]), match[2]), (int(match[3]), match[4])


def _prompt_annotation_for_range(mol_id, chain, first, last, refresh_function=None, selection_cleanup=None):
  try:
    if mol_id not in model_molecule_list():
      raise ValueError("The annotation model is no longer available.")
    snapshot = _model_coordinate_snapshot(mol_id)
    start, end = (chain, *first), (chain, *last)
    specs = [s for s in snapshot if s[0] == chain]
    if start not in specs or end not in specs:
      raise ValueError("Both endpoints must exist in the chosen chain.")
    if not _residue_is_polymer(mol_id, *start) or not _residue_is_polymer(mol_id, *end):
      raise ValueError("Region notes require polymer endpoints in the same chain.")
    lower, upper = sorted((specs.index(start), specs.index(end)))
    start, end = specs[lower], specs[upper]
    def target_for(spec):
      return {"mol_id": mol_id, "chain_id": spec[0], "resno": spec[1], "ins_code": spec[2],
              "residue_name": snapshot[spec]["name"]}
    target = target_for(start)
    target["end"] = target_for(end)
    label = _annotation_entry_label(_new_annotation_entry(target, "", "", "", 1))
    generic_title_author_note_entry("Add region note", label, "", _annotation_default_author(), "",
      "Add region note", lambda title, author, note: _append_annotation_for_target(
        target, title, author, note, refresh_function), close_function=selection_cleanup)
    return 1
  except (ValueError, RuntimeError) as error:
    if callable(selection_cleanup):
      selection_cleanup()
    info_dialog(str(error))
    return 0


def _annotation_click_spec_for_model(spec, mol_id):
  """Use the annotation model to repair GTK4's unset molecule field."""
  if not isinstance(spec, (list, tuple)) or len(spec) not in (6, 7):
    raise ValueError("Coot returned an unrecognised atom pick.")
  index = 1 if len(spec) == 7 else 0
  clicked_model = _coot_int_arg(spec[index])
  if clicked_model != -1:
    if clicked_model != mol_id or mol_id not in model_molecule_list():
      raise ValueError("Click an endpoint in the annotation model.")
    return list(spec)
  chain, number, ins = _click_spec_chain_id(spec), _click_spec_res_no(spec), _click_spec_ins_code(spec)
  atom_name, alt = _click_spec_atom_name(spec).strip(), _click_spec_alt_conf(spec).strip()
  if not chain or number is False or not atom_name:
    raise ValueError("Coot returned an incomplete atom pick.")
  if mol_id not in model_molecule_list():
    raise ValueError("The annotation model is no longer available.")
  displayed = getattr(coot, "mol_is_displayed", None)
  if callable(displayed) and not displayed(mol_id):
    raise ValueError("Display the annotation model before picking the range.")
  if not residue_exists_qm(mol_id, chain, number, ins):
    raise ValueError("The clicked atom is not in the displayed annotation model.")
  atoms = coot.residue_info_py(mol_id, chain, number, ins) or []
  if not any(a[0][0].strip() == atom_name and str(a[0][1] or "").strip() == alt for a in atoms):
    raise ValueError("The clicked atom is not in the displayed annotation model.")
  repaired = list(spec)
  repaired[index] = mol_id
  return repaired


def _clear_annotation_range_pick():
  state = dict(RESIDUE_ANNOTATION_RANGE_PICK_STATE)
  RESIDUE_ANNOTATION_RANGE_PICK_STATE.clear()
  if state.get("picking"):
    clear_picks = getattr(coot, "clear_pending_picks", None)
    if callable(clear_picks):
      try:
        clear_picks()
      except Exception as error:
        print("Could not cancel annotation endpoint picking:", error)
  for obj in state.get("markers", []):
    try:
      coot.close_generic_object(obj)
    except Exception as error:
      print("Could not clear annotation endpoint marker:", error)
  window = state.get("window")
  if window is not None:
    window.destroy()
  draw = getattr(coot, "graphics_draw", None)
  if state.get("markers") and callable(draw):
    draw()


def _mark_annotation_endpoint(spec, colour):
  required = ("new_generic_object_number", "to_generic_object_add_point",
              "set_display_generic_object", "close_generic_object")
  if any(not callable(getattr(coot, name, None)) for name in required):
    return
  atoms = coot.residue_info_py(_click_spec_imol(spec), _click_spec_chain_id(spec),
                             _click_spec_res_no(spec), _click_spec_ins_code(spec)) or []
  candidates = [a for a in atoms if a[0][0].strip() == "CA"]
  if not candidates:
    candidates = [a for a in atoms if a[0][0].strip() == "P"] or atoms
  if not candidates:
    return
  atom = max(candidates, key=lambda a: (not bool(a[0][1]), float(a[1][0])))
  obj = coot.new_generic_object_number("Annotation range endpoint")
  RESIDUE_ANNOTATION_RANGE_PICK_STATE.setdefault("markers", []).append(obj)
  try:
    coot.to_generic_object_add_point(obj, colour, 10, *atom[2])
    coot.set_display_generic_object(obj, 1)
    draw = getattr(coot, "graphics_draw", None)
    if callable(draw):
      draw()
  except Exception as error:
    print("Could not draw annotation endpoint marker:", error)
    coot.close_generic_object(obj)
    RESIDUE_ANNOTATION_RANGE_PICK_STATE["markers"].remove(obj)


def _start_annotation_range_pick(mol_id, refresh_function=None):
  cancel_script = globals().get("_cancel_script_clicks")
  if callable(cancel_script):
    cancel_script()
  _clear_annotation_range_pick()
  token = object()
  RESIDUE_ANNOTATION_RANGE_PICK_STATE.update(token=token, picking=True, markers=[])
  def current():
    return RESIDUE_ANNOTATION_RANGE_PICK_STATE.get("token") is token
  def cleanup():
    if current():
      _clear_annotation_range_pick()
  def decode(spec):
    spec = _annotation_click_spec_for_model(spec, mol_id)
    if not _residue_is_polymer(mol_id, _click_spec_chain_id(spec),
                               _click_spec_res_no(spec), _click_spec_ins_code(spec)):
      raise ValueError("Region endpoints must be polymer residues.")
    return spec
  def show_click(spec, suffix=""):
    add_status_bar_text(_navigation_status_bar_label(mol_id, _click_spec_chain_id(spec),
      _click_spec_res_no(spec), _click_spec_ins_code(spec)).replace("Current residue:", "Clicked residue:", 1) + suffix)
  def reject(error, callback):
    print("Annotation range picking:", error)
    add_status_bar_text(str(error))
    if current():
      coot.user_defined_click_py(1, callback)
  def last_clicked(first, raw_last):
    if not current():
      return
    try:
      first = decode(first)
      last = decode(raw_last)
      chain = _click_spec_chain_id(first)
      if _click_spec_chain_id(last) != chain:
        raise ValueError("Click the end residue in the same chain as the start residue.")
    except ValueError as error:
      reject(error, lambda last: last_clicked(first, last))
      return
    show_click(last)
    _mark_annotation_endpoint(last, "orange")
    RESIDUE_ANNOTATION_RANGE_PICK_STATE["picking"] = False
    window = RESIDUE_ANNOTATION_RANGE_PICK_STATE.pop("window", None)
    if window is not None:
      window.destroy()
    try:
      _prompt_annotation_for_range(mol_id, chain, (_click_spec_res_no(first), _click_spec_ins_code(first)),
                                  (_click_spec_res_no(last), _click_spec_ins_code(last)), refresh_function,
                                  selection_cleanup=cleanup)
    except Exception:
      cleanup()
      raise
  def first_clicked(raw_first):
    if not current():
      return
    try:
      first = decode(raw_first)
    except ValueError as error:
      reject(error, first_clicked)
      return
    show_click(first, " | Click the end residue (orange). Tab cancels.")
    _mark_annotation_endpoint(first, "cyan")
    coot.user_defined_click_py(1, lambda last: last_clicked(first, last))
  add_status_bar_text("Region note: click start (cyan), then end (orange). Tab cancels.")
  try:
    coot.user_defined_click_py(1, first_clicked)
  except Exception:
    cleanup()
    raise


def prompt_add_annotation_for_region(mol_id=None, refresh_function=None):
  if mol_id is None:
    mol_id = _annotation_active_molecule_or_status()
  if mol_id is None:
    return 0
  def enter_range():
    residue = active_residue()
    chain = residue[1] if residue and residue[0] == mol_id else "A"
    def submitted(chain_text, range_text):
      try:
        first, last = _annotation_parse_range(range_text)
      except ValueError as error:
        info_dialog(str(error))
        return 0
      return _prompt_annotation_for_range(mol_id, chain_text.strip(), first, last, refresh_function)
    generic_double_entry("Add region note", "Chain", chain, "Residue range", "75-80",
                         "Continue", submitted)
  action_button_dialog("Add region note", [["Click two endpoints", lambda: _start_annotation_range_pick(mol_id, refresh_function)],
                                           ["Enter chain / range...", enter_range]])
  return 1


def _annotation_note_choice_label(entry, entry_index, max_length=48):
  author_text = _annotation_text_or_empty(entry.get("author")) or "Unknown author"
  modified_text = _annotation_text_or_empty(entry.get("modified_utc")) or "Unknown time"
  preview_text = _annotation_entry_display_title(entry, fallback_from_note=True, max_length=max_length)
  if preview_text:
    return f"{entry_index}. {author_text} - {modified_text} - {preview_text}"
  return f"{entry_index}. {author_text} - {modified_text}"


def _update_annotation_entry(annotation_entry, title_text, author_text, note_text, refresh_function=None):
  global RESIDUE_ANNOTATION_LAST_AUTHOR
  note_fields = _annotation_prepare_note_fields(title_text, author_text, note_text)
  if note_fields is None:
    return 0
  normalized_title, normalized_author, normalized_note = note_fields
  RESIDUE_ANNOTATION_LAST_AUTHOR = normalized_author
  annotation_entry["title"] = normalized_title
  annotation_entry["author"] = normalized_author
  annotation_entry["note"] = normalized_note
  annotation_entry["modified_utc"] = _annotation_now_utc()
  if callable(refresh_function):
    refresh_function()
  add_status_bar_text("Updated residue annotation")
  return 1


def _delete_annotation_entry(mol_id, annotation_entry, refresh_function=None):
  entries = _annotation_entries_for_molecule(mol_id)
  try:
    entries.remove(annotation_entry)
  except ValueError:
    info_dialog("That residue annotation is no longer present.")
    return 0
  if callable(refresh_function):
    refresh_function()
  add_status_bar_text("Deleted residue annotation")
  return 1


def _prompt_edit_annotation_entry(group, annotation_entry, refresh_function=None):
  generic_title_author_note_entry(
    "Edit residue annotation",
    group["dialog_label"],
    _annotation_text_or_empty(annotation_entry.get("title")),
    _annotation_text_or_empty(annotation_entry.get("author")) or _annotation_default_author(),
    _annotation_text_or_empty(annotation_entry.get("note")),
    "Save changes",
    lambda title_text, author_text, note_text: _update_annotation_entry(
      annotation_entry,
      title_text,
      author_text,
      note_text,
      refresh_function,
    ),
  )
  return 1


def _choose_annotation_entry_action(group, action_name, action_function):
  entries = list(group.get("entries") or [])
  if not entries:
    info_dialog("No note is selected for this residue.")
    return 0
  if len(entries) == 1:
    return action_function(entries[0])
  button_list = []
  for entry_index, annotation_entry in enumerate(entries, start=1):
    button_list.append(
      (
        _annotation_note_choice_label(annotation_entry, entry_index),
        lambda current_entry=annotation_entry: action_function(current_entry),
      )
    )
  action_button_dialog(f"{action_name} - {group['dialog_label']}", button_list)
  return 1


def _annotation_jump_to_group(mol_id, group):
  entry = group["entries"][0]
  if entry.get("target_kind") == "region":
    specs = _annotation_region_specs(mol_id, entry)
    if not specs:
      _clear_odd_helical_feature_highlights(RESIDUE_ANNOTATION_REGION_HIGHLIGHTS)
      info_dialog("The annotated region endpoints are missing or ambiguous in this model.")
      return 0
    _clear_odd_helical_feature_highlights()
    if not _highlight_odd_helical_range({"mol_id": mol_id, "residue_specs": specs},
                                        RESIDUE_ANNOTATION_REGION_HIGHLIGHTS):
      return 0
    chain_id, resno, ins_code = specs[len(specs)//2]
    result = _go_to_navigation_residue(mol_id, chain_id, resno, ins_code,
                                      _residue_serial_number(mol_id, chain_id, resno, ins_code))
    add_status_bar_text(f"{group['dialog_label']} - {_annotation_entry_display_title(entry)}")
    return result
  _clear_odd_helical_feature_highlights(RESIDUE_ANNOTATION_REGION_HIGHLIGHTS)
  chain_id, resno, ins_code, residue_name_here = group["key"][:4]
  residue_point = residue_centre_py(mol_id, chain_id, resno, ins_code)
  if not isinstance(residue_point, (list, tuple)) or len(residue_point) != 3:
    residue_point = None
  jump_entry = [group["dialog_label"]]
  if residue_point is not None:
    jump_entry.extend(residue_point)
  else:
    jump_entry.extend([0.0, 0.0, 0.0])

  if _residue_is_polymer(mol_id, chain_id, resno, ins_code):
    jump_entry.append(
      {
        "type": "polymer_residue",
        "mol_id": mol_id,
        "chain_id": chain_id,
        "resno": resno,
        "ins_code": ins_code or "",
        "serial_number": _residue_serial_number(mol_id, chain_id, resno, ins_code),
      }
    )
  elif residue_point is None:
    add_status_bar_text(group["dialog_label"])
    return 0

  _activate_interesting_entry(jump_entry)
  return 1


def _annotation_groups_near_point(mol_id, point_xyz, radius, groups=None, positions=None):
  nearby_groups = []
  try:
    radius_sq = float(radius) * float(radius)
  except Exception:
    return nearby_groups
  groups = _annotation_groups_for_molecule(mol_id) if groups is None else groups
  all_specs = None
  if any(g.get("entries") and g["entries"][0].get("target_kind") == "region" for g in groups):
    all_specs = list(positions) if positions is not None else _all_residue_specs_for_colouring(mol_id)
  for group in groups:
    entry = group.get("entries", [{}])[0]
    if entry.get("target_kind") == "region":
      specs = group.get("region_specs")
      if specs is None:
        specs = _annotation_region_specs(mol_id, entry, all_specs)
    else:
      specs = [group["key"][:3]]
    distances = []
    for spec in specs:
      point = residue_centre_py(mol_id, *spec) if positions is None else positions.get(tuple(spec))
      if isinstance(point, (list, tuple)) and len(point) == 3:
        distances.append(_distance_sq(point_xyz, point))
    if not distances:
      continue
    distance_sq = min(distances)
    if distance_sq > radius_sq:
      continue
    nearby_group = dict(group)
    nearby_group["distance"] = math.sqrt(distance_sq)
    nearby_groups.append(nearby_group)
  nearby_groups.sort(key=lambda group: (group.get("distance", float("inf")), _annotation_group_sort_key(group)))
  return nearby_groups


def _annotation_groups_signature(groups):
  return tuple(
    (
      group["key"],
      len(group.get("entries") or []),
      _annotation_text_or_empty(group["entries"][-1].get("modified_utc")) if group.get("entries") else "",
      round(float(group.get("distance", -1.0)), 2),
    )
    for group in groups
  )


def _annotation_browser_title(mol_id):
  source_display_name = _annotation_source_display_name(mol_id)
  if source_display_name:
    return (
      f"Residue annotations for molecule #{mol_id}\n"
      f"Source: {source_display_name}; use the widget to export an annotated mmCIF"
    )
  return f"Residue annotations for molecule #{mol_id}\n(No accessible source file path)"


def _annotation_nearby_browser_title(mol_id, radius):
  source_display_name = _annotation_source_display_name(mol_id)
  radius_text = f"{float(radius):.1f} A"
  if source_display_name:
    return (
      f"Nearby residue annotations within {radius_text} of the rotation centre for molecule #{mol_id}\n"
      f"Source: {source_display_name}"
    )
  return (
    f"Nearby residue annotations within {radius_text} of the rotation centre for molecule #{mol_id}\n"
    "(No accessible source file path)"
  )


def _annotation_active_residue_header_text(mol_id):
  residue = active_residue()
  if not residue:
    return f"Active residue: none in molecule #{mol_id}"
  active_mol_id, chain_id, resno, ins_code = residue[0:4]
  residue_name_here = residue_name(active_mol_id, chain_id, resno, ins_code) or ""
  residue_label = _odd_residue_dialog_label(chain_id, resno, ins_code, residue_name_here)
  if active_mol_id == mol_id:
    return f"Active residue: {residue_label}"
  return f"Active residue: {residue_label} (molecule #{active_mol_id})"


def residue_annotations_gui():
  initial_mol_id = _annotation_active_molecule_or_status()
  if initial_mol_id is None:
    return 0

  _ensure_residue_annotations_loaded(initial_mol_id)

  if "Gtk" not in globals():
    browser_data = _annotation_browser_data(initial_mol_id)
    grouped = browser_data["groups"]
    labels = [_annotation_group_preview_text(group) for group in grouped]
    message = browser_data["title"]
    if labels:
      message += "\n\n" + "\n".join(labels)
    info_dialog(message)
    return len(grouped)

  window = Gtk.Window()
  window.set_title("Coot")
  window.set_default_size(640, 620)

  outer_vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=8)
  action_row_top = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
  action_row_bottom = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
  navigation_row = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
  title_label = Gtk.Label(label=_annotation_browser_title(initial_mol_id))
  list_scrolled = Gtk.ScrolledWindow()
  list_box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6)
  active_residue_label = Gtk.Label(label=_annotation_active_residue_header_text(initial_mol_id))
  detail_header = Gtk.Label(label="Selected notes: none")
  detail_scrolled = Gtk.ScrolledWindow()
  detail_box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=4)
  add_button = Gtk.Button(label="Add note for active residue")
  add_region_button = Gtk.Button(label="Add region note...")
  edit_button = Gtk.Button(label="Edit selected note...")
  delete_button = Gtk.Button(label="Delete selected note...")
  export_button = Gtk.Button(label="Export annotated mmCIF...")
  export_markdown_button = Gtk.Button(label="Export Markdown table...")
  clear_button = Gtk.Button(label="Clear all notes...")
  refresh_button = Gtk.Button(label="Refresh list")
  prev_button = Gtk.Button(label="<--Prev")
  next_button = Gtk.Button(label="Next-->")
  close_button = Gtk.Button(label="Close")

  title_label.set_margin_start(12)
  title_label.set_margin_end(12)
  title_label.set_margin_top(12)
  title_label.set_margin_bottom(4)
  title_label.set_wrap(True)
  title_label.set_xalign(0.0)
  list_scrolled.set_margin_start(12)
  list_scrolled.set_margin_end(12)
  list_scrolled.set_margin_bottom(8)
  list_scrolled.set_hexpand(True)
  list_scrolled.set_vexpand(True)
  try:
    list_scrolled.set_has_frame(True)
  except Exception:
    pass
  list_box.set_valign(Gtk.Align.START)
  list_scrolled.set_child(list_box)
  active_residue_label.set_margin_start(12)
  active_residue_label.set_margin_end(12)
  active_residue_label.set_margin_bottom(2)
  active_residue_label.set_xalign(0.0)
  detail_header.set_margin_start(12)
  detail_header.set_margin_end(12)
  detail_header.set_margin_bottom(2)
  detail_header.set_xalign(0.0)
  detail_scrolled.set_margin_start(12)
  detail_scrolled.set_margin_end(12)
  detail_scrolled.set_margin_bottom(8)
  detail_scrolled.set_hexpand(True)
  detail_scrolled.set_vexpand(True)
  try:
    detail_scrolled.set_has_frame(True)
  except Exception:
    pass
  detail_box.set_valign(Gtk.Align.START)
  detail_scrolled.set_child(detail_box)
  action_row_top.set_halign(Gtk.Align.CENTER)
  action_row_top.set_margin_start(12)
  action_row_top.set_margin_end(12)
  action_row_top.set_margin_bottom(4)
  action_row_bottom.set_halign(Gtk.Align.CENTER)
  action_row_bottom.set_margin_start(12)
  action_row_bottom.set_margin_end(12)
  action_row_bottom.set_margin_bottom(4)
  navigation_row.set_halign(Gtk.Align.CENTER)
  navigation_row.set_margin_start(12)
  navigation_row.set_margin_end(12)
  navigation_row.set_margin_bottom(4)
  close_button.set_margin_start(12)
  close_button.set_margin_end(12)
  close_button.set_margin_bottom(12)
  close_button.set_halign(Gtk.Align.CENTER)

  glib_module = None

  try:
    glib_module = _coot_gui_repository_module("GLib")
  except Exception:
    glib_module = None

  class AnnotationBrowserState:
    def __init__(self, mol_id):
      self.current_mol_id = mol_id
      self.grouped_annotations = []
      self.group_index_by_key = {}
      self.selected_group_key = None
      self.current_index = None
      self.selected_button = None
      self.row_buttons = []
      self.active_residue_timer_id = None
      self.last_active_residue_header = None

    def mol_id(self):
      return self.current_mol_id

    def stop_polling(self):
      _clear_annotation_range_pick()
      _clear_odd_helical_feature_highlights(RESIDUE_ANNOTATION_REGION_HIGHLIGHTS)
      if glib_module is not None and self.active_residue_timer_id is not None:
        try:
          glib_module.source_remove(self.active_residue_timer_id)
        except Exception:
          pass
        self.active_residue_timer_id = None
      return None

    def close_window(self, *_args):
      self.stop_polling()
      window.destroy()
      return False

    def start_polling(self):
      if glib_module is not None:
        try:
          self.active_residue_timer_id = glib_module.timeout_add(
            RESIDUE_ANNOTATION_POLL_INTERVAL_MS,
            self.poll_active_residue_header,
          )
        except Exception:
          self.active_residue_timer_id = None
      return None

    def selected_group(self):
      if self.selected_group_key is None:
        return None
      group_index = self.group_index_by_key.get(self.selected_group_key)
      if group_index is None or group_index >= len(self.grouped_annotations):
        return None
      return self.grouped_annotations[group_index]

    def current_mol_id_or_dialog(self):
      mol_id = self.mol_id()
      if mol_id in model_molecule_list():
        return mol_id
      info_dialog("The current annotation molecule is no longer available.")
      return None

    def selected_group_or_dialog(self):
      selected_group = self.selected_group()
      if selected_group is None:
        info_dialog("Select a residue annotation first.")
      return selected_group

    def set_selected_button(self, button):
      if self.selected_button is not None:
        try:
          self.selected_button.remove_css_class("suggested-action")
        except Exception:
          pass
      self.selected_button = button
      if button is not None:
        try:
          button.add_css_class("suggested-action")
        except Exception:
          pass

    def clear_list_box_children(self):
      _annotation_clear_box_children(list_box)

    def update_selection_button_state(self):
      has_selection = self.current_index is not None
      edit_button.set_sensitive(has_selection)
      delete_button.set_sensitive(has_selection)

    def set_active_residue_header_text(self, header_text):
      active_residue_label.set_text(header_text)
      self.last_active_residue_header = header_text

    def update_active_residue_header(self):
      header_text = _annotation_active_residue_header_text(self.mol_id())
      if header_text != self.last_active_residue_header:
        self.set_active_residue_header_text(header_text)
      return None

    def clear_selected_group(self):
      _clear_odd_helical_feature_highlights(RESIDUE_ANNOTATION_REGION_HIGHLIGHTS)
      self.current_index = None
      self.selected_group_key = None
      detail_header.set_text("Selected notes: none")
      _annotation_set_detail_placeholder(detail_box, "Select a residue annotation thread to view its notes.")
      self.set_selected_button(None)
      self.update_selection_button_state()

    def show_group(self, entry_index, jump_to_residue=True):
      if not self.grouped_annotations:
        return None
      entry_index = entry_index % len(self.grouped_annotations)
      group = self.grouped_annotations[entry_index]
      self.current_index = entry_index
      self.selected_group_key = group["key"]
      detail_header.set_text(f"Selected notes: {group['dialog_label']}")
      _annotation_render_group_detail_widgets(detail_box, group)
      button = self.row_buttons[entry_index] if entry_index < len(self.row_buttons) else None
      self.set_selected_button(button)
      self.update_selection_button_state()
      if jump_to_residue:
        _annotation_jump_to_group(self.mol_id(), group)
      return group

    def activate_prev(self, *_args):
      if not self.grouped_annotations:
        return None
      if self.current_index is None:
        return self.show_group(len(self.grouped_annotations) - 1)
      return self.show_group(self.current_index - 1)

    def activate_next(self, *_args):
      if not self.grouped_annotations:
        return None
      if self.current_index is None:
        return self.show_group(0)
      return self.show_group(self.current_index + 1)

    def render_groups(self, preferred_key=None, jump_to_selection=False):
      browser_data = _annotation_browser_data(self.mol_id())
      self.grouped_annotations[:] = browser_data["groups"]
      self.group_index_by_key = {}
      title_label.set_text(browser_data["title"])
      self.clear_selected_group()
      self.row_buttons.clear()
      self.clear_list_box_children()
      self.set_active_residue_header_text(_annotation_active_residue_header_text(self.mol_id()))

      if not self.grouped_annotations:
        empty_label = Gtk.Label(label="No residue annotations for this molecule.")
        empty_label.set_xalign(0.0)
        list_box.append(empty_label)
        prev_button.set_sensitive(False)
        next_button.set_sensitive(False)
        return None

      selected_index = None
      for entry_index, group in enumerate(self.grouped_annotations):
        self.group_index_by_key[group["key"]] = entry_index
        button = Gtk.Button(label=_annotation_group_preview_text(group))
        button.connect("clicked", lambda _button, idx=entry_index: self.show_group(idx))
        self.row_buttons.append(button)
        list_box.append(button)
        if preferred_key is not None and group["key"] == preferred_key:
          selected_index = entry_index

      prev_button.set_sensitive(True)
      next_button.set_sensitive(True)
      if selected_index is not None:
        self.show_group(selected_index, jump_to_residue=jump_to_selection)
      else:
        self.update_active_residue_header()
      return None

    def refresh_groups(self, preferred_key=None):
      self.render_groups(
        preferred_key if preferred_key is not None else self.selected_group_key,
        jump_to_selection=False,
      )
      return None

    def poll_active_residue_header(self):
      model_molecules = model_molecule_list()
      active_mol_id = _annotation_active_molecule()
      target_mol_id = active_mol_id
      if target_mol_id is None and self.mol_id() in model_molecules:
        target_mol_id = self.mol_id()
      if target_mol_id is not None and target_mol_id != self.mol_id():
        _clear_annotation_range_pick()
        self.current_mol_id = target_mol_id
        _ensure_residue_annotations_loaded(target_mol_id)
        self.render_groups()
      self.update_active_residue_header()
      return True

    def add_note(self, *_args):
      mol_id = self.current_mol_id_or_dialog()
      if mol_id is None:
        return None
      return prompt_add_annotation_for_active_residue(mol_id, self.refresh_groups)

    def add_region_note(self, *_args):
      mol_id = self.current_mol_id_or_dialog()
      if mol_id is not None:
        return prompt_add_annotation_for_region(mol_id, self.refresh_groups)

    def edit_selected_note(self, *_args):
      if self.current_mol_id_or_dialog() is None:
        return None
      selected_group = self.selected_group_or_dialog()
      if selected_group is None:
        return None
      return _choose_annotation_entry_action(
        selected_group,
        "Edit note",
        lambda annotation_entry: _prompt_edit_annotation_entry(
          selected_group,
          annotation_entry,
          lambda preferred_key=selected_group["key"]: self.refresh_groups(preferred_key),
        ),
      )

    def delete_selected_note(self, *_args):
      selected_mol_id = self.current_mol_id_or_dialog()
      if selected_mol_id is None:
        return None
      selected_group = self.selected_group_or_dialog()
      if selected_group is None:
        return None
      return _choose_annotation_entry_action(
        selected_group,
        "Delete note",
        lambda annotation_entry: generic_confirm_dialog(
          "Delete residue annotation",
          "Delete the selected note from this residue thread?",
          "Cancel",
          lambda: None,
          "Delete note",
          lambda: _delete_annotation_entry(
            selected_mol_id,
            annotation_entry,
            lambda preferred_key=selected_group["key"]: self.refresh_groups(preferred_key),
          ),
        ),
      )

    def export_annotated_mmcif(self, *_args):
      mol_id = self.current_mol_id_or_dialog()
      if mol_id is None:
        return None
      default_output_path = _default_annotation_export_path(mol_id)
      generic_single_entry_with_file_browse(
        "Export annotated mmCIF",
        default_output_path,
        "Choose annotated mmCIF output path",
        "Export annotated mmCIF",
        lambda output_path: _export_residue_annotations_with_confirmation(mol_id, output_path),
      )
      return None

    def export_markdown_table(self, *_args):
      mol_id = self.current_mol_id_or_dialog()
      if mol_id is None:
        return None
      default_output_path = _default_annotation_markdown_export_path(mol_id)
      generic_single_entry_with_file_browse(
        "Export annotation Markdown",
        default_output_path,
        "Choose annotation Markdown output path",
        "Export Markdown table",
        lambda output_path: _export_residue_annotations_markdown_with_confirmation(mol_id, output_path),
      )
      return None

    def clear_all_notes_for_molecule(self, target_mol_id):
      _annotation_entries_for_molecule(target_mol_id).clear()
      if target_mol_id == self.mol_id():
        self.refresh_groups()
      add_status_bar_text("Cleared all residue annotations for the current molecule")
      return None

    def clear_all_notes(self, *_args):
      selected_mol_id = self.current_mol_id_or_dialog()
      if selected_mol_id is None:
        return None
      if not _annotation_entries_for_molecule(selected_mol_id):
        info_dialog("No residue annotations are currently stored for this molecule.")
        return None
      generic_confirm_dialog(
        "Clear all residue annotations",
        "Delete all residue annotations for the current molecule?\n"
        "This only affects the in-memory/browser state until you export an annotated mmCIF.",
        "Cancel",
        lambda: None,
        "Clear all notes",
        lambda: self.clear_all_notes_for_molecule(selected_mol_id),
      )
      return None

  browser = AnnotationBrowserState(initial_mol_id)

  add_button.connect("clicked", browser.add_note)
  add_region_button.connect("clicked", browser.add_region_note)
  edit_button.connect("clicked", browser.edit_selected_note)
  delete_button.connect("clicked", browser.delete_selected_note)
  export_button.connect("clicked", browser.export_annotated_mmcif)
  export_markdown_button.connect("clicked", browser.export_markdown_table)
  clear_button.connect("clicked", browser.clear_all_notes)
  refresh_button.connect("clicked", lambda *_args: browser.refresh_groups())
  prev_button.connect("clicked", browser.activate_prev)
  next_button.connect("clicked", browser.activate_next)
  close_button.connect("clicked", browser.close_window)
  window.connect("close-request", lambda *_args: browser.stop_polling() or False)
  export_button.set_sensitive(True)

  browser.render_groups()
  browser.start_polling()

  action_row_top.append(add_button)
  action_row_top.append(add_region_button)
  action_row_top.append(edit_button)
  action_row_top.append(delete_button)
  action_row_bottom.append(export_button)
  action_row_bottom.append(export_markdown_button)
  action_row_bottom.append(clear_button)
  action_row_bottom.append(refresh_button)
  navigation_row.append(prev_button)
  navigation_row.append(next_button)

  outer_vbox.append(title_label)
  outer_vbox.append(action_row_top)
  outer_vbox.append(action_row_bottom)
  outer_vbox.append(list_scrolled)
  outer_vbox.append(active_residue_label)
  outer_vbox.append(detail_header)
  outer_vbox.append(detail_scrolled)
  outer_vbox.append(navigation_row)
  outer_vbox.append(close_button)
  window.set_child(outer_vbox)
  window.present()
  return len(browser.grouped_annotations)


def nearby_residue_annotations_gui(radius=RESIDUE_ANNOTATION_NEARBY_RADIUS):
  initial_mol_id = _annotation_active_molecule_or_status()
  if initial_mol_id is None:
    return 0

  try:
    nearby_radius = float(radius)
  except Exception:
    nearby_radius = float(RESIDUE_ANNOTATION_NEARBY_RADIUS)
  if nearby_radius <= 0.0:
    nearby_radius = float(RESIDUE_ANNOTATION_NEARBY_RADIUS)

  _ensure_residue_annotations_loaded(initial_mol_id)

  if "Gtk" not in globals():
    groups = _annotation_groups_near_point(initial_mol_id, _rotation_centre_xyz(), nearby_radius)
    labels = [_annotation_group_nearby_preview_text(group) for group in groups]
    message = _annotation_nearby_browser_title(initial_mol_id, nearby_radius)
    if labels:
      message += "\n\n" + "\n".join(labels)
    else:
      message += "\n\nNo residue annotations are currently near the rotation centre."
    info_dialog(message)
    return len(groups)

  window = Gtk.Window()
  window.set_title("Coot")
  window.set_default_size(620, 560)

  outer_vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=8)
  action_row = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
  navigation_row = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
  title_label = Gtk.Label(label=_annotation_nearby_browser_title(initial_mol_id, nearby_radius))
  context_label = Gtk.Label(label="")
  list_scrolled = Gtk.ScrolledWindow()
  list_box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6)
  detail_header = Gtk.Label(label="Nearby notes: none")
  detail_scrolled = Gtk.ScrolledWindow()
  detail_box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=4)
  refresh_button = Gtk.Button(label="Refresh nearby notes")
  prev_button = Gtk.Button(label="<--Prev")
  next_button = Gtk.Button(label="Next-->")
  close_button = Gtk.Button(label="Close")

  title_label.set_margin_start(12)
  title_label.set_margin_end(12)
  title_label.set_margin_top(12)
  title_label.set_margin_bottom(4)
  title_label.set_wrap(True)
  title_label.set_xalign(0.0)
  context_label.set_margin_start(12)
  context_label.set_margin_end(12)
  context_label.set_margin_bottom(4)
  context_label.set_wrap(True)
  context_label.set_xalign(0.0)
  list_scrolled.set_margin_start(12)
  list_scrolled.set_margin_end(12)
  list_scrolled.set_margin_bottom(8)
  list_scrolled.set_hexpand(True)
  list_scrolled.set_vexpand(True)
  try:
    list_scrolled.set_has_frame(True)
  except Exception:
    pass
  list_box.set_valign(Gtk.Align.START)
  list_scrolled.set_child(list_box)
  detail_header.set_margin_start(12)
  detail_header.set_margin_end(12)
  detail_header.set_margin_bottom(2)
  detail_header.set_xalign(0.0)
  detail_scrolled.set_margin_start(12)
  detail_scrolled.set_margin_end(12)
  detail_scrolled.set_margin_bottom(8)
  detail_scrolled.set_hexpand(True)
  detail_scrolled.set_vexpand(True)
  try:
    detail_scrolled.set_has_frame(True)
  except Exception:
    pass
  detail_box.set_valign(Gtk.Align.START)
  detail_scrolled.set_child(detail_box)
  action_row.set_halign(Gtk.Align.CENTER)
  action_row.set_margin_start(12)
  action_row.set_margin_end(12)
  action_row.set_margin_bottom(4)
  navigation_row.set_halign(Gtk.Align.CENTER)
  navigation_row.set_margin_start(12)
  navigation_row.set_margin_end(12)
  navigation_row.set_margin_bottom(4)
  close_button.set_margin_start(12)
  close_button.set_margin_end(12)
  close_button.set_margin_bottom(12)
  close_button.set_halign(Gtk.Align.CENTER)

  glib_module = None
  try:
    glib_module = _coot_gui_repository_module("GLib")
  except Exception:
    glib_module = None

  class NearbyAnnotationBrowserState:
    def __init__(self, mol_id, nearby_radius):
      self.current_mol_id = mol_id
      self.radius = float(nearby_radius)
      self.grouped_annotations = []
      self.group_index_by_key = {}
      self.selected_group_key = None
      self.current_index = None
      self.selected_button = None
      self.row_buttons = []
      self.timer_id = None
      self.last_signature = None
      self.last_rotation_centre = None
      self.cached_annotation_entries = None
      self.cached_annotation_groups = []
      self.last_nearby_query = None
      self.cached_positions = {}
      self.positions_updated_at = -float("inf")
      self.positions_epoch = None

    def mol_id(self):
      return self.current_mol_id

    def stop_polling(self):
      if glib_module is not None and self.timer_id is not None:
        try:
          glib_module.source_remove(self.timer_id)
        except Exception:
          pass
        self.timer_id = None
      return None

    def close_window(self, *_args):
      self.stop_polling()
      window.destroy()
      return False

    def start_polling(self):
      if glib_module is not None:
        try:
          self.timer_id = glib_module.timeout_add(
            RESIDUE_ANNOTATION_POLL_INTERVAL_MS,
            self.poll_nearby_groups,
          )
        except Exception:
          self.timer_id = None
      return None

    def clear_list_box_children(self):
      _annotation_clear_box_children(list_box)

    def set_selected_button(self, button):
      if self.selected_button is not None:
        try:
          self.selected_button.remove_css_class("suggested-action")
        except Exception:
          pass
      self.selected_button = button
      if button is not None:
        try:
          button.add_css_class("suggested-action")
        except Exception:
          pass

    def selected_group(self):
      if self.selected_group_key is None:
        return None
      group_index = self.group_index_by_key.get(self.selected_group_key)
      if group_index is None or group_index >= len(self.grouped_annotations):
        return None
      return self.grouped_annotations[group_index]

    def clear_selected_group(self):
      self.current_index = None
      self.selected_group_key = None
      detail_header.set_text("Nearby notes: none")
      _annotation_set_detail_placeholder(detail_box, "No nearby residue annotation thread is selected.")
      self.set_selected_button(None)

    def update_context_text(self):
      count = len(self.grouped_annotations)
      if count == 0:
        context_label.set_text(
          f"No residue annotation threads are currently within {self.radius:.1f} A of the rotation centre."
        )
      else:
        context_label.set_text(
          f"{count} residue annotation {'thread' if count == 1 else 'threads'} within "
          f"{self.radius:.1f} A of the rotation centre."
        )
      return None

    def show_group(self, entry_index):
      if not self.grouped_annotations:
        return None
      entry_index = entry_index % len(self.grouped_annotations)
      group = self.grouped_annotations[entry_index]
      self.current_index = entry_index
      self.selected_group_key = group["key"]
      detail_header.set_text(f"Nearby notes: {group['dialog_label']}")
      _annotation_render_group_detail_widgets(detail_box, group)
      button = self.row_buttons[entry_index] if entry_index < len(self.row_buttons) else None
      self.set_selected_button(button)
      return group

    def activate_prev(self, *_args):
      if not self.grouped_annotations:
        return None
      if self.current_index is None:
        return self.show_group(len(self.grouped_annotations) - 1)
      return self.show_group(self.current_index - 1)

    def activate_next(self, *_args):
      if not self.grouped_annotations:
        return None
      if self.current_index is None:
        return self.show_group(0)
      return self.show_group(self.current_index + 1)

    def render_groups(self, groups, preferred_key=None):
      self.grouped_annotations[:] = list(groups or [])
      self.group_index_by_key = {}
      title_label.set_text(_annotation_nearby_browser_title(self.mol_id(), self.radius))
      self.clear_selected_group()
      self.row_buttons.clear()
      self.clear_list_box_children()

      if not self.grouped_annotations:
        empty_label = Gtk.Label(label="No residue annotations are currently near the rotation centre.")
        empty_label.set_xalign(0.0)
        list_box.append(empty_label)
        prev_button.set_sensitive(False)
        next_button.set_sensitive(False)
        self.update_context_text()
        return None

      selected_index = None
      for entry_index, group in enumerate(self.grouped_annotations):
        self.group_index_by_key[group["key"]] = entry_index
        button = Gtk.Button(label=_annotation_group_nearby_preview_text(group))
        button.connect("clicked", lambda _button, idx=entry_index: self.show_group(idx))
        self.row_buttons.append(button)
        list_box.append(button)
        if preferred_key is not None and group["key"] == preferred_key:
          selected_index = entry_index

      prev_button.set_sensitive(True)
      next_button.set_sensitive(True)
      self.update_context_text()
      if selected_index is None:
        selected_index = 0
      self.show_group(selected_index)
      return None

    def refresh_nearby_groups(self, force=False):
      model_molecules = model_molecule_list()
      active_mol_id = _annotation_active_molecule()
      target_mol_id = active_mol_id
      if target_mol_id is None and self.mol_id() in model_molecules:
        target_mol_id = self.mol_id()
      if target_mol_id is None:
        self.grouped_annotations[:] = []
        self.last_signature = None
        self.last_rotation_centre = None
        self.last_nearby_query = None
        self.clear_selected_group()
        self.clear_list_box_children()
        empty_label = Gtk.Label(label="No active model molecule is available for nearby notes.")
        empty_label.set_xalign(0.0)
        list_box.append(empty_label)
        prev_button.set_sensitive(False)
        next_button.set_sensitive(False)
        context_label.set_text("No active model molecule is available for nearby residue annotations.")
        return None

      if target_mol_id != self.mol_id():
        self.current_mol_id = target_mol_id
        _ensure_residue_annotations_loaded(target_mol_id)
        force = True

      rotation_centre = _rotation_centre_xyz()
      entries = _annotation_entries_for_molecule(target_mol_id)
      if force or entries != self.cached_annotation_entries:
        self.cached_annotation_entries = [dict(entry) for entry in entries]
        self.cached_annotation_groups = _annotation_groups_for_molecule(target_mol_id)
        force = True
      now = time.monotonic()
      epoch = _MODEL_EDIT_EPOCHS.get(target_mol_id, 0)
      # Script edits invalidate immediately. Native edits/Undo outside this
      # script are picked up by the bounded fallback refresh.
      if force or epoch != self.positions_epoch or now - self.positions_updated_at >= 1.5:
        positions = {}
        has_regions = any(group.get("entries") and group["entries"][0].get("target_kind") == "region"
                          for group in self.cached_annotation_groups)
        bulk_centres = (_coordinate_snapshot_centres(_model_coordinate_snapshot(target_mol_id))
                        if len(self.cached_annotation_groups) >= 64 or has_regions else None)
        if has_regions:
          needed = set()
          for group in self.cached_annotation_groups:
            entry = group.get("entries", [{}])[0]
            if entry.get("target_kind") == "region":
              group["region_specs"] = _annotation_region_specs(target_mol_id, entry, bulk_centres)
              needed.update(group["region_specs"])
            else:
              needed.add(group["key"][:3])
          positions.update((spec, point) for spec, point in bulk_centres.items() if spec in needed)
        for group in self.cached_annotation_groups:
          key = group["key"][:3]
          if key not in positions:
            position = bulk_centres.get(key) if bulk_centres is not None else residue_centre_py(target_mol_id, *key)
            positions[key] = tuple(position) if isinstance(position, (list, tuple)) and len(position) == 3 else None
        self.cached_positions = positions
        self.positions_updated_at = now
        self.positions_epoch = epoch
      positions = self.cached_positions
      query = (target_mol_id, tuple(rotation_centre), self.radius, tuple(positions.items()))
      if not force and query == self.last_nearby_query:
        return None
      self.last_nearby_query = query
      groups = _annotation_groups_near_point(target_mol_id, rotation_centre, self.radius,
                                             groups=self.cached_annotation_groups, positions=positions)
      signature = _annotation_groups_signature(groups)
      moved = (
        self.last_rotation_centre is None
        or _distance_sq(rotation_centre, self.last_rotation_centre) > 0.25 * 0.25
      )
      if force or moved or signature != self.last_signature:
        self.last_signature = signature
        self.last_rotation_centre = list(rotation_centre)
        self.render_groups(groups, preferred_key=self.selected_group_key)
      else:
        self.grouped_annotations[:] = list(groups or [])
        self.update_context_text()
      return None

    def poll_nearby_groups(self):
      self.refresh_nearby_groups(force=False)
      return True

  browser = NearbyAnnotationBrowserState(initial_mol_id, nearby_radius)

  refresh_button.connect("clicked", lambda *_args: browser.refresh_nearby_groups(force=True))
  prev_button.connect("clicked", browser.activate_prev)
  next_button.connect("clicked", browser.activate_next)
  close_button.connect("clicked", browser.close_window)
  window.connect("close-request", lambda *_args: browser.stop_polling() or False)

  action_row.append(refresh_button)
  navigation_row.append(prev_button)
  navigation_row.append(next_button)

  outer_vbox.append(title_label)
  outer_vbox.append(context_label)
  outer_vbox.append(action_row)
  outer_vbox.append(list_scrolled)
  outer_vbox.append(detail_header)
  outer_vbox.append(detail_scrolled)
  outer_vbox.append(navigation_row)
  outer_vbox.append(close_button)
  window.set_child(outer_vbox)

  browser.refresh_nearby_groups(force=True)
  browser.start_polling()
  window.present()
  return len(browser.grouped_annotations)


def _generate_smart_local_extra_restraints_for_mol(
  mol_id,
  distance_cutoff=3.7,
  show_start_message=True,
):
  if mol_id not in model_molecule_list():
    add_status_bar_text("No active model")
    return None
  if show_start_message:
    add_status_bar_text("Generating smart local restraints...")
    graphics_draw()
  try:
    distance_cutoff = float(distance_cutoff)
  except (TypeError, ValueError):
    distance_cutoff = 3.7
  if not math.isfinite(distance_cutoff) or distance_cutoff <= 0:
    info_dialog("Restraint distance cutoff must be finite and positive.")
    return None
  max_sequence_separation = 10
  restraint_esd = 0.05
  backbone_donor_names = {"N"}
  backbone_acceptor_names = {"O", "OXT"}

  sequence_index_by_residue = {}
  next_sequence_index_by_chain = {}
  ordered_residue_specs_by_chain = {}
  model_snapshot = _model_coordinate_snapshot(mol_id)
  for residue_spec_list in model_snapshot:
    residue_spec = _coot_residue_spec_from_spec(residue_spec_list)
    if residue_spec is None:
      continue
    chain_id, resno, ins_code = residue_spec
    sequence_index_by_residue[residue_spec] = next_sequence_index_by_chain.get(chain_id, 0)
    next_sequence_index_by_chain[chain_id] = next_sequence_index_by_chain.get(chain_id, 0) + 1
    ordered_residue_specs_by_chain.setdefault(chain_id, []).append(residue_spec)

  def atom_distance(atom_1, atom_2):
    dx = atom_1["position"][0] - atom_2["position"][0]
    dy = atom_1["position"][1] - atom_2["position"][1]
    dz = atom_1["position"][2] - atom_2["position"][2]
    return math.sqrt(dx*dx + dy*dy + dz*dz)

  def residue_atoms(residue_spec):
    return model_snapshot[residue_spec]["atoms"]

  def parsed_atom_record(atom):
    position = residue_atom_to_position(atom)
    if not position:
      return None
    return {
      "name": residue_atom_to_atom_name(atom),
      "alt_conf": residue_atom2alt_conf(atom),
      "position": position,
    }

  def residue_record(residue_spec):
    parsed_atoms = []
    backbone_hbond_atoms = []
    for atom in residue_atoms(residue_spec):
      parsed_atom = parsed_atom_record(atom)
      if not parsed_atom:
        continue
      parsed_atoms.append(parsed_atom)
      atom_name = parsed_atom["name"].strip()
      if atom_name in backbone_donor_names or atom_name in backbone_acceptor_names:
        backbone_hbond_atoms.append(parsed_atom)
    if not parsed_atoms:
      return None
    is_polymer = model_snapshot[residue_spec]["name"] in POLYMER_RESIDUE_NAMES
    return {
      "spec": residue_spec,
      "atoms": parsed_atoms,
      "backbone_hbond_atoms": backbone_hbond_atoms,
      "polymer": is_polymer,
    }

  def sequence_separation(residue_spec_1, residue_spec_2):
    if residue_spec_1[0] != residue_spec_2[0]:
      return None
    index_1 = sequence_index_by_residue.get(residue_spec_1)
    index_2 = sequence_index_by_residue.get(residue_spec_2)
    if index_1 is None or index_2 is None:
      return None
    return abs(index_1 - index_2)

  def is_backbone_hbond_candidate(atom_name_1, atom_name_2):
    atom_name_1 = atom_name_1.strip()
    atom_name_2 = atom_name_2.strip()
    return (
      (atom_name_1 in backbone_donor_names and atom_name_2 in backbone_acceptor_names)
      or (atom_name_1 in backbone_acceptor_names and atom_name_2 in backbone_donor_names)
    )

  def ordered_pair_key(spec_1, spec_2):
    if spec_2 < spec_1:
      return (spec_2, spec_1)
    return (spec_1, spec_2)

  def add_restraints_between_records(record_1, record_2, backbone_only=False):
    nonlocal added_restraints
    if not record_1 or not record_2:
      return
    atoms_1 = record_1["backbone_hbond_atoms"] if backbone_only else record_1["atoms"]
    atoms_2 = record_2["backbone_hbond_atoms"] if backbone_only else record_2["atoms"]
    if backbone_only and not (record_1["polymer"] and record_2["polymer"]):
      return
    for atom_1 in atoms_1:
      atom_name_1 = atom_1["name"]
      alt_conf_1 = atom_1["alt_conf"]
      for atom_2 in atoms_2:
        atom_name_2 = atom_2["name"]
        alt_conf_2 = atom_2["alt_conf"]
        if backbone_only and not is_backbone_hbond_candidate(atom_name_1, atom_name_2):
          continue
        if alt_conf_1 and alt_conf_2 and alt_conf_1 != alt_conf_2:
          continue
        restraint_key = tuple(sorted((
          (*record_1["spec"], atom_name_1, alt_conf_1),
          (*record_2["spec"], atom_name_2, alt_conf_2))))
        if restraint_key in seen_restraints:
          continue
        distance = atom_distance(atom_1, atom_2)
        if distance is None or distance >= distance_cutoff:
          continue
        pending_restraints.append((
          mol_id,
          record_1["spec"][0], record_1["spec"][1], record_1["spec"][2], atom_name_1, alt_conf_1,
          record_2["spec"][0], record_2["spec"][1], record_2["spec"][2], atom_name_2, alt_conf_2,
          distance, restraint_esd
        ))
        seen_restraints.add(restraint_key)
        added_restraints += 1

  added_restraints = 0
  pending_restraints = []
  seen_restraints = set()
  processed_backbone_hbond_pairs = set()
  residue_records = {}
  for chain_id, residue_specs in ordered_residue_specs_by_chain.items():
    for residue_spec in residue_specs:
      residue_records[residue_spec] = residue_record(residue_spec)

  for chain_id, residue_specs in ordered_residue_specs_by_chain.items():
    for index, residue_spec_1 in enumerate(residue_specs):
      record_1 = residue_records.get(residue_spec_1)
      if not record_1:
        continue
      upper_index = min(len(residue_specs), index + max_sequence_separation + 1)
      for neighbour_index in range(index + 1, upper_index):
        residue_spec_2 = residue_specs[neighbour_index]
        add_restraints_between_records(record_1, residue_records.get(residue_spec_2))

  hbond_atoms = [(record["spec"], atom) for record in residue_records.values()
                 if record and record["polymer"] for atom in record["backbone_hbond_atoms"]]
  hbond_index = _spatial_index(hbond_atoms, lambda item: item[1]["position"], distance_cutoff)
  all_residue_specs = list(residue_records.keys())
  for residue_spec_1 in all_residue_specs:
    record_1 = residue_records.get(residue_spec_1)
    if not record_1 or not record_1["polymer"]:
      continue
    nearby_residue_specs = {spec for atom in record_1["backbone_hbond_atoms"]
                            for spec, _other in _spatial_query(hbond_index, atom["position"], distance_cutoff)}
    for nearby_spec in nearby_residue_specs:
      residue_spec_2 = _coot_residue_spec_from_spec(nearby_spec)
      if residue_spec_2 is None:
        continue
      if residue_spec_1 == residue_spec_2:
        continue
      pair_key = ordered_pair_key(residue_spec_1, residue_spec_2)
      if pair_key in processed_backbone_hbond_pairs:
        continue
      processed_backbone_hbond_pairs.add(pair_key)
      separation = sequence_separation(residue_spec_1, residue_spec_2)
      if separation is not None and separation <= max_sequence_separation:
        continue
      add_restraints_between_records(record_1, residue_records.get(residue_spec_2), backbone_only=True)

  # Finish all discovery before replacing existing restraints.
  delete_all_extra_restraints(mol_id)
  for args in pending_restraints:
    add_extra_geman_mcclure_restraint(*args)
  set_show_extra_restraints(mol_id, 0)
  set_show_extra_restraints(mol_id, 1)
  set_show_extra_distance_restraints(1)
  summary = "Smart local extra restraints generated: {added} added".format(
      added=added_restraints
    )
  add_status_bar_text(summary)
  print(summary)
  return added_restraints


def generate_smart_local_extra_restraints():
  residue = _active_residue_or_status()
  if not residue:
    return None
  return _generate_smart_local_extra_restraints_for_mol(residue[0], show_start_message=True)


def generate_smart_local_extra_restraints_with_cutoff(distance_cutoff):
  residue = _active_residue_or_status()
  if not residue:
    return None
  try:
    parsed_cutoff = float(distance_cutoff)
  except ValueError:
    info_dialog("Minimum interatomic distance must be a number")
    return None
  if parsed_cutoff <= 0.0:
    info_dialog("Minimum interatomic distance must be greater than 0")
    return None
  return _generate_smart_local_extra_restraints_for_mol(
    residue[0], distance_cutoff=parsed_cutoff, show_start_message=True
  )


def prompt_generate_smart_local_extra_restraints():
  """Prompt for the local distance cutoff before building smart restraints."""
  generic_single_entry(
    "Minimum interatomic distance for smart restraints (A)",
    "3.7",
    "Generate smart self restraints",
    generate_smart_local_extra_restraints_with_cutoff,
  )


def flip_active_peptide():
  # The installed coot_fitting helper may call the removed coot.closest_atom_raw()
  # name; use the supported Python bindings directly instead.
  active_atom = coot.closest_atom_simple_py() if hasattr(coot, "closest_atom_simple_py") else False
  if not active_atom:
    add_status_bar_text("No active residue")
    return None
  imol = active_atom[0]
  atom_spec = coot.closest_atom_raw_py() if hasattr(coot, "closest_atom_raw_py") else False
  if isinstance(atom_spec, (list, tuple)) and atom_spec and isinstance(atom_spec[0], (list, tuple)):
    atom_spec = atom_spec[0]
  if not isinstance(atom_spec, (list, tuple)) or len(atom_spec) < 6:
    add_status_bar_text("Could not identify the peptide at the screen centre")
    return None
  chain_id = atom_spec[1]
  res_no = atom_spec[2]
  ins_code = atom_spec[3]
  atom_name = atom_spec[4]
  alt_conf = atom_spec[5]
  if str(atom_name).strip() == "N":
    res_no = res_no - 1
  return coot.pepflip(imol, chain_id, res_no, ins_code, alt_conf)


def add_water_and_refine():
  mol_id = place_water_in_active_molecule()
  if mol_id is None:
    return None
  placed_residue = _active_residue_or_status()
  if not placed_residue:
    return None
  placed_water_atom = _residue_first_atom_for_refinement_guard(
    placed_residue[0],
    placed_residue[1],
    placed_residue[2],
    placed_residue[3],
  )
  if imol_refinement_map()==-1:
    add_status_bar_text("Water added; no refinement map set")
    return None
  previous_immediate_replacement = refinement_immediate_replacement_state()
  try:
    set_refinement_immediate_replacement(1)
    if hasattr(coot, "rsr_sphere_refine"):
      rsr_sphere_refine()
      accept_result = accept_moving_atoms_py()
      if accept_result is False:
        add_status_bar_text("Water added, but local refinement could not be accepted")
        return None
    else:
      residue = _active_residue_or_status()
      if not residue:
        return None
      ch_id = residue[1]
      resno = residue[2]
      altloc = residue[5]
      refine_zone(mol_id, ch_id, resno, resno, altloc)
      accept_regularizement()
  finally:
    set_refinement_immediate_replacement(previous_immediate_replacement)
  _restore_water_if_refined_too_far(
    placed_residue,
    placed_water_atom,
    WATER_REFINE_MAX_DISPLACEMENT,
  )
  return mol_id


def place_water_in_active_molecule():
  residue = _active_residue_or_status()
  if not residue:
    return None
  mol_id = residue[0]
  try:
    set_pointer_atom_molecule(mol_id)
  except Exception:
    pass
  try:
    set_go_to_atom_molecule(mol_id)
  except Exception:
    pass
  place_typed_atom_at_pointer("Water")
  return mol_id


def _residue_first_atom_for_refinement_guard(mol_id, ch_id, resno, ins_code):
  """Return the first atom in a residue with enough metadata to restore xyz."""
  atom_info = residue_info_py(mol_id, ch_id, resno, ins_code)
  if not isinstance(atom_info, list):
    return None
  for atom in atom_info:
    try:
      atom_name = atom[0][0]
      alt_conf = atom[0][1]
      xyz = list(atom[2])
    except Exception:
      continue
    if atom_name and len(xyz) == 3:
      return {
        "atom_name": atom_name,
        "alt_conf": alt_conf or "",
        "xyz": xyz,
      }
  return None


def _restore_atom_xyz(mol_id, ch_id, resno, ins_code, atom_name, alt_conf, xyz):
  """Move one atom back to a stored Cartesian position."""
  for attribute_name, value in zip(("x", "y", "z"), xyz):
    set_atom_attribute(mol_id, ch_id, resno, ins_code, atom_name, alt_conf, attribute_name, value)


def _restore_water_if_refined_too_far(placed_residue, placed_water_atom, max_displacement):
  """Undo pathological water refinement jumps by restoring the original placement."""
  if not placed_residue or not placed_water_atom:
    return None
  refined_water_atom = _residue_first_atom_for_refinement_guard(
    placed_residue[0],
    placed_residue[1],
    placed_residue[2],
    placed_residue[3],
  )
  if not refined_water_atom:
    return None
  if _distance_sq(placed_water_atom["xyz"], refined_water_atom["xyz"]) <= (max_displacement * max_displacement):
    return None
  _restore_atom_xyz(
    placed_residue[0],
    placed_residue[1],
    placed_residue[2],
    placed_residue[3],
    refined_water_atom["atom_name"],
    refined_water_atom["alt_conf"],
    placed_water_atom["xyz"],
  )
  add_status_bar_text(
    "Water moved more than {0:.1f} A during refinement; restored original placement".format(
      max_displacement,
    )
  )
  return 1


def undo_symmetry_view_safe():
  if hasattr(coot, "undo_symmetry_view"):
    return coot.undo_symmetry_view()
  add_status_bar_text("Undo symmetry view is unavailable in this build")
  return None


def narrow_clipping_symmetric():
  increase_clipping_front()
  decrease_clipping_back()


def widen_clipping_symmetric():
  decrease_clipping_front()
  increase_clipping_back()


def _map_radius_controls(map_id, allow_inferred=False):
  native_type = getattr(coot, "is_EM_map", None)
  if callable(native_type):
    is_em = bool(native_type(map_id))
  else:
    if not allow_inferred:
      # Older APIs cannot report user-overridden/PANDDA map classification.
      # Keep both radius settings synchronized rather than update the wrong one.
      return None
    # Match Coot's default P1/orthogonal-cell classification, not filename.
    sg = str(coot.space_group_py(map_id) or "").replace(" ", "").upper()
    parameters = map_cell(map_id)
    if not parameters or len(parameters) < 6 or not sg:
      return None
    is_em = sg == "P1" and all(abs(math.radians(float(angle)-90.0)) < 0.0001 for angle in parameters[3:6])
  return (get_map_radius_em, set_map_radius_em) if is_em else (get_map_radius, set_map_radius)


def _current_map_radius(map_id=None):
  if map_id is None:
    map_id = _scrollable_map_or_status()
  controls = _map_radius_controls(map_id, allow_inferred=True) if map_id is not None else None
  return controls[0]() if controls else get_map_radius()


def _set_map_radius_both(radius, map_id=None):
  """Compatibility name: change only the radius used by the target map."""
  if map_id is None:
    map_id = _scrollable_map_or_status()
  controls = _map_radius_controls(map_id) if map_id is not None else None
  for getter, setter in ([controls] if controls else
                         [(get_map_radius, set_map_radius), (get_map_radius_em, set_map_radius_em)]):
    if abs(float(getter()) - float(radius)) > 1e-5:
      setter(radius)


def increase_map_radius():
  current_radius = _current_map_radius()
  _set_map_radius_both(current_radius + 2.0)


def decrease_map_radius():
  current_radius = _current_map_radius()
  new_radius = max(2.0, current_radius - 2.0)
  _set_map_radius_both(new_radius)


def _proportional_editing_radius_supported():
  return all(callable(getattr(coot, name, None)) for name in (
    "increase_proportional_editing_radius", "decrease_proportional_editing_radius",
    "refinement_already_ongoing_p"))


def _show_proportional_editing_radius_unavailable():
  info_dialog(
    "This Coot build does not expose proportional editing radius controls yet."
  )
  return None


def _proportional_editing_radius_step_size():
  # Upstream 44e657b (revision-count 19684) changed 1.0 to 0.3 A. HEAD
  # builds report "1.x-pre", so check their revision count before the version.
  revision_function = getattr(coot, "git_revision_count", None)
  if callable(revision_function) and int(revision_function()) >= 19684:
    return 0.3
  version_function = getattr(coot, "coot_version", None)
  version_string = str(version_function()) if callable(version_function) else ""
  version = re.search(r"(\d+)\.(\d+)(?:\.(\d+))?", version_string)
  return 0.3 if version and tuple(int(part or 0) for part in version.groups()) >= (1, 3, 0) else 1.0


def _apply_pending_proportional_editing_radius(from_timer=True):
  """Coot ignores radius steps until refinement restraints have been created."""
  global PROPORTIONAL_EDITING_RADIUS
  state = PROPORTIONAL_EDITING_RADIUS_STATE
  target = state["pending"]
  if target is None:
    if from_timer:
      state["source"] = None
    return False
  try:
    if not coot.refinement_already_ongoing_p():
      return True
    # There is no native getter. Track only steps made while refinement exists;
    # Coot's initial radius is zero, not the previous widget's assumed 1 A.
    step = _proportional_editing_radius_step_size()
    while state["applied"] < target - step / 2.0:
      coot.increase_proportional_editing_radius()
      state["applied"] = round(state["applied"] + step, 6)
    while state["applied"] > target + step / 2.0:
      coot.decrease_proportional_editing_radius()
      state["applied"] = max(0.0, round(state["applied"] - step, 6))
    state["pending"] = None
    add_status_bar_text("Proportional editing radius: %.1f A" % target)
  except Exception as error:
    state["pending"] = None
    PROPORTIONAL_EDITING_RADIUS = state["applied"]
    info_dialog("Unable to change the proportional editing radius.\n\n%r" % (error,))
  if from_timer:
    state["source"] = None
  return False


def _set_proportional_editing_radius_value(requested_radius):
  global PROPORTIONAL_EDITING_RADIUS
  if not _proportional_editing_radius_supported():
    return _show_proportional_editing_radius_unavailable()
  try:
    requested_radius = float(requested_radius)
  except (TypeError, ValueError):
    info_dialog("Radius must be a number")
    return False
  if not math.isfinite(requested_radius) or not 0.0 <= requested_radius <= 1000.0:
    info_dialog("Proportional editing radius must be between 0 and 1000 A")
    return False
  step = _proportional_editing_radius_step_size()
  target = round(round(requested_radius / step) * step, 6)
  state = PROPORTIONAL_EDITING_RADIUS_STATE
  previous = PROPORTIONAL_EDITING_RADIUS
  PROPORTIONAL_EDITING_RADIUS = target
  if target == state["applied"]:
    state["pending"] = None
    add_status_bar_text("Proportional editing radius: %.1f A" % target)
    return True
  state["pending"] = target
  try:
    if coot.refinement_already_ongoing_p():
      _apply_pending_proportional_editing_radius(from_timer=False)
      return state["applied"] == target
    if state["source"] is None:
      glib = _coot_gui_repository_module("GLib")
      if glib is None:
        raise RuntimeError("GLib is needed to apply this setting when refinement starts")
      state["source"] = glib.timeout_add(250, _apply_pending_proportional_editing_radius)
  except Exception as error:
    state["pending"] = None
    PROPORTIONAL_EDITING_RADIUS = previous
    info_dialog("Could not schedule the proportional editing radius.\n\n%r" % (error,))
    return False
  add_status_bar_text("Proportional editing radius: %.1f A (applies when refinement starts)" % target)
  return True


def _apply_proportional_editing_radius_step(direction):
  step = _proportional_editing_radius_step_size()
  target = min(1000.0, max(0.0, PROPORTIONAL_EDITING_RADIUS + (step if direction >= 0 else -step)))
  if _set_proportional_editing_radius_value(target):
    return target
  return None


def increase_proportional_editing_radius_with_status():
  return _apply_proportional_editing_radius_step(+1.0)


def decrease_proportional_editing_radius_with_status():
  return _apply_proportional_editing_radius_step(-1.0)


def set_proportional_editing_radius():
  if not _proportional_editing_radius_supported():
    return _show_proportional_editing_radius_unavailable()
  current_radius_string = "%.1f" % float(PROPORTIONAL_EDITING_RADIUS)
  step = _proportional_editing_radius_step_size()
  def set_proportional_editing_radius_from_text(text):
    return int(bool(_set_proportional_editing_radius_value(text)))
  generic_single_entry(
    "Proportional editing radius (A; %.1f A steps; below 2 A = off)" % step,
    current_radius_string,
    "Set proportional editing radius",
    set_proportional_editing_radius_from_text,
  )


def active_map_surface_displayed():
  map_id = _scrollable_map_or_status()
  if map_id is None:
    return False
  return (
    map_is_displayed(map_id)
    and MAP_SURFACE_DISPLAY_STATE.get(map_id, False)
  )


def _active_surface_map_or_status():
  """Return the active map if it is currently shown as a solid surface."""
  if not active_map_surface_displayed():
    add_status_bar_text("Active map is not in solid-surface mode")
    return None
  return scroll_wheel_map()


def _adjust_active_map_surface_opacity(delta):
  """Adjust the opacity of the active solid-surface map."""
  map_id = _active_surface_map_or_status()
  if map_id is None:
    return None
  current_opacity = MAP_SURFACE_OPACITY_STATE.get(
    map_id, get_solid_density_surface_opacity(map_id)
  )
  new_opacity = min(1.0, max(0.0, current_opacity + delta))
  MAP_SURFACE_OPACITY_STATE[map_id] = new_opacity
  set_solid_density_surface_opacity(map_id, new_opacity)


def increase_active_map_surface_opacity():
  """Increase active-map solid-surface opacity in 0.05 steps."""
  return _adjust_active_map_surface_opacity(0.05)


def decrease_active_map_surface_opacity():
  """Decrease active-map solid-surface opacity in 0.05 steps."""
  return _adjust_active_map_surface_opacity(-0.05)


def display_only_active_map():
  """Show only the active map, cycling if one map is already displayed alone."""
  map_ids = map_molecule_list()
  if not map_ids:
    add_status_bar_text("No maps loaded")
    return None
  active_map = scroll_wheel_map()
  if active_map not in map_ids:
    for map_id in map_molecule_list():
      if (map_is_displayed(map_id)==1) and (map_id!=active_map):
        set_scroll_wheel_map(map_id)
        set_scrollable_map(map_id)
      else:
        set_scroll_wheel_map(map_ids[0])
        set_scrollable_map(map_ids[0])
    active_map=scroll_wheel_map()
  if not map_is_displayed(active_map):
    set_map_displayed(active_map,1)
  displayed_maps_count=0
  for map_id in map_molecule_list():
    displayed_maps_count=displayed_maps_count+map_is_displayed(map_id)
    if (map_is_displayed(map_id)==1) and (map_id!=active_map):
      set_map_displayed(map_id,0)
    if map_is_displayed(map_id):
      displayed_map=map_id
  if displayed_maps_count==1:
    index_displayed=map_molecule_list().index(active_map)
    try:
      next_map=map_molecule_list()[index_displayed+1]
    except IndexError:
      next_map=map_molecule_list()[0]
    set_map_displayed(active_map,0)
    set_map_displayed(next_map,1)
    set_scroll_wheel_map(next_map)
    set_scrollable_map(next_map)
  for map_id in map_molecule_list():
    if map_is_displayed(map_id):
      set_scrollable_map(map_id)
      set_scroll_wheel_map(map_id) #New

def hide_active_mol():
  residue = _active_residue_or_status()
  if not residue:
    return None
  mol_id=residue[0]
  set_mol_displayed(mol_id,0)

def display_only_active():
  mol_id_active = _active_model_molecule_for_map_tools("display-only toggle")
  if mol_id_active is None:
    return None
  displayed_mols_count=0
  for mol_id in model_molecule_list():
    displayed_mols_count=displayed_mols_count+mol_is_displayed(mol_id)
    if (mol_is_displayed(mol_id)==1) and (mol_id!=mol_id_active):
      set_mol_displayed(mol_id,0)
    elif (mol_is_displayed(mol_id)==0) and (mol_id==mol_id_active):
      set_mol_displayed(mol_id,1)
    if mol_is_displayed(mol_id):
      displayed_mol=mol_id
  if displayed_mols_count==1:
    index_displayed=model_molecule_list().index(mol_id_active)
    try: 
      next_mol=model_molecule_list()[index_displayed+1]
    except IndexError:
      next_mol=model_molecule_list()[0]
    set_mol_displayed(displayed_mol,0)
    set_mol_displayed(next_mol,1)


def _map_level_summary(map_id):
  """Return the current map level in absolute and sigma units."""
  contour_absolute = get_contour_level_absolute(map_id)
  contour_sigma = get_contour_level_in_sigma(map_id)
  return contour_absolute, contour_sigma


def _post_map_level_status(map_id, prefix="set to"):
  """Report the current map level in a consistent status-bar format."""
  contour_absolute, contour_sigma = _map_level_summary(map_id)
  add_status_bar_text(
    "Map #{0} {1} level {2:.3f}, {3:.1f} sigma".format(
      map_id,
      prefix,
      contour_absolute,
      contour_sigma,
    )
  )


def _tq_status(message):
  print(message)
  add_status_bar_text(message)


def _tq_residue_info(imol, chain_id, resno, ins_code):
  chain_id, resno, ins_code = _coot_residue_spec_tuple(chain_id, resno, ins_code)
  if resno is None:
    return []
  return residue_info_py(imol, chain_id, resno, ins_code) or []


def _tq_residue_name(imol, chain_id, resno, ins_code):
  chain_id, resno, ins_code = _coot_residue_spec_tuple(chain_id, resno, ins_code)
  if resno is None:
    return ""
  return residue_name(imol, chain_id, resno, ins_code) or ""


def _tq_residue_centre(imol, chain_id, resno, ins_code):
  chain_id, resno, ins_code = _coot_residue_spec_tuple(chain_id, resno, ins_code)
  if resno is None:
    return None
  centre = residue_centre_py(imol, chain_id, resno, ins_code)
  if isinstance(centre, (list, tuple)) and len(centre) == 3:
    return [float(centre[0]), float(centre[1]), float(centre[2])]
  return None


def _tq_map_id_or_default(map_id=None):
  if map_id is not None:
    return map_id
  return _scrollable_map_or_status()


def _tq_map_sigma(map_id):
  try:
    sigma_value = float(map_sigma_py(map_id))
    if math.isfinite(sigma_value) and sigma_value > 0.0:
      return sigma_value
  except Exception:
    pass
  try:
    contour_absolute = float(get_contour_level_absolute(map_id))
    contour_sigma = float(get_contour_level_in_sigma(map_id))
    if math.isfinite(contour_absolute) and math.isfinite(contour_sigma) and abs(contour_sigma) > 1.0e-9:
      sigma_value = abs(contour_absolute / contour_sigma)
      if sigma_value > 0.0:
        return sigma_value
  except Exception:
    pass
  return None


def _tq_valid_model_molecule(mol_id):
  checker = globals().get("valid_model_molecule_qm")
  if callable(checker):
    try:
      return checker(mol_id) == 1 or checker(mol_id) is True
    except Exception:
      pass
  return mol_id in model_molecule_list()


def _tq_valid_map_molecule(map_id):
  checker = globals().get("valid_map_molecule_qm")
  if callable(checker):
    try:
      return checker(map_id) == 1 or checker(map_id) is True
    except Exception:
      pass
  return map_id in map_molecule_list()


def _tq_chain_ids(mol_id):
  getter = globals().get("chain_ids")
  if callable(getter):
    try:
      return getter(mol_id) or []
    except Exception:
      pass
  return []


def _tq_post_map_level_status(map_id, prefix="set to"):
  _post_map_level_status(map_id, prefix)


def _tq_active_residue_target():
  residue = _active_residue_or_status()
  if not residue or len(residue) < 4:
    return None
  parsed = _coot_residue_spec_from_spec(residue[1:4])
  if parsed is None:
    _tq_status("Thresholded Q-residue: could not parse the active residue")
    return None
  return (residue[0], parsed[0], parsed[1], parsed[2])


def _tq_is_heavy_atom(atom_name, element):
  element = _coot_text_arg(element).strip().upper()
  atom_name = _coot_text_arg(atom_name).strip().upper()
  if element:
    return element not in _THRESHOLDED_Q_HYDROGEN_ELEMENTS
  return not atom_name.startswith("H")


def _tq_residue_label(chain_id, resno, ins_code, residue_name_here=""):
  residue_id = f"{chain_id}:{resno}{ins_code or ''}"
  if residue_name_here:
    return f"{residue_id} {_display_residue_name(residue_name_here)}"
  return residue_id


def _tq_gaussian_enclosed_mass_fraction(radius_over_sigma):
  if radius_over_sigma <= 0.0:
    return 0.0
  scaled_radius = radius_over_sigma / math.sqrt(2.0)
  return (
    math.erf(scaled_radius)
    - math.sqrt(2.0 / math.pi) * radius_over_sigma * math.exp(-(radius_over_sigma ** 2) / 2.0)
  )


@lru_cache(maxsize=64)
def _tq_support_radius(sigma, enclosed_mass_percent):
  target_fraction = max(0.0, min(0.999999, float(enclosed_mass_percent) / 100.0))
  if target_fraction <= 0.0:
    return 0.0
  if target_fraction >= 0.999999:
    return 5.0 * float(sigma)
  low = 0.0
  high = 8.0
  for _ in range(60):
    mid = 0.5 * (low + high)
    if _tq_gaussian_enclosed_mass_fraction(mid) < target_fraction:
      low = mid
    else:
      high = mid
  return 0.5 * (low + high) * float(sigma)


def _tq_residue_atom_records(mol_id, residue_spec, heavy_atoms_only=True, snapshot=None):
  parsed = _coot_residue_spec_from_spec(residue_spec)
  if parsed is None:
    return []
  chain_id, resno, ins_code = parsed
  cached = snapshot.get(parsed) if snapshot is not None else None
  residue_name_here = cached["name"] if cached is not None else _tq_residue_name(mol_id, chain_id, resno, ins_code)
  atom_records = []
  atoms = cached["atoms"] if cached is not None else _tq_residue_info(mol_id, chain_id, resno, ins_code)
  for atom in atoms:
    try:
      atom_name = _coot_text_arg(atom[0][0]).strip()
      alt_conf = _coot_text_arg(atom[0][1]).strip()
      occupancy = float(atom[1][0])
      b_factor = atom[1][1]
      element = _coot_text_arg(atom[1][2]).strip().upper()
      xyz = [float(atom[2][0]), float(atom[2][1]), float(atom[2][2])]
    except Exception:
      continue
    if not atom_name:
      continue
    if not math.isfinite(occupancy) or not all(math.isfinite(value) for value in xyz):
      raise RuntimeError(f"Invalid atom coordinates or occupancy at {chain_id}:{resno}{ins_code}.")
    if heavy_atoms_only and not _tq_is_heavy_atom(atom_name, element):
      continue
    atom_records.append(
      {
        "residue_spec": parsed,
        "residue_name": residue_name_here,
        "name": atom_name,
        "alt_conf": alt_conf,
        "atom_key": (parsed, atom_name, alt_conf),
        "occupancy": occupancy,
        "b_factor": b_factor,
        "element": element,
        "xyz": xyz,
      }
    )
  return atom_records


def _tq_score_atom_records(atom_records, backbone_only=False):
  if not backbone_only:
    return list(atom_records)
  return [atom for atom in atom_records if atom["name"] in _THRESHOLDED_Q_BACKBONE_ATOM_NAMES]


def _tq_chain_residue_specs(mol_id, chain_id):
  residue_specs = []
  seen_specs = set()
  try:
    n_residues = chain_n_residues(chain_id, mol_id)
  except Exception:
    return residue_specs
  for serial_number in range(n_residues):
    try:
      resno = seqnum_from_serial_number(mol_id, chain_id, serial_number)
      ins_code = _coot_text_arg(insertion_code_from_serial_number(mol_id, chain_id, serial_number))
    except Exception:
      continue
    residue_spec = _coot_residue_spec_tuple(chain_id, resno, ins_code)
    if residue_spec[1] is None or residue_spec in seen_specs:
      continue
    residue_specs.append(residue_spec)
    seen_specs.add(residue_spec)
  return residue_specs


def _tq_residue_specs_near_residue_centres(
    mol_id,
    seed_residue_specs,
    radius=THRESHOLDED_Q_LOCAL_RADIUS,
    max_residues=THRESHOLDED_Q_LOCAL_MAX_RESIDUES,
    snapshot=None,
):
  snapshot = _model_coordinate_snapshot(mol_id) if snapshot is None else snapshot
  centres = _coordinate_snapshot_centres(snapshot)
  seed_specs = []
  seed_centres = []
  for residue_spec in seed_residue_specs:
    parsed = _coot_residue_spec_from_spec(residue_spec)
    if parsed is None:
      continue
    seed_specs.append(parsed)
    centre = centres.get(parsed)
    if centre is not None:
      seed_centres.append(centre)
  if not seed_specs:
    return []
  if not seed_centres:
    return list(seed_specs)

  radius_sq = float(radius) * float(radius)
  nearby_entries = []
  seen_specs = set(seed_specs)
  for seed_spec in seed_specs:
    nearby_entries.append((0.0, seed_spec))

  centre_index = _spatial_index(centres.items(), lambda item: item[1], max(float(radius), 1.0))
  for seed_centre in seed_centres:
    for residue_spec, centre in _spatial_query(centre_index, seed_centre, float(radius)):
      if residue_spec not in seen_specs:
        nearest_seed_distance_sq = min(_distance_sq(centre, point) for point in seed_centres)
        nearby_entries.append((nearest_seed_distance_sq, residue_spec))
        seen_specs.add(residue_spec)

  nearby_entries.sort(key=lambda item: (item[0], item[1][0], item[1][1], item[1][2]))
  return [entry[1] for entry in (nearby_entries if max_residues is None else nearby_entries[:int(max_residues)])]


def _tq_context_atom_records(mol_id, score_residue_specs, context_radius, snapshot=None,
                             influence_radius=0.0):
  snapshot = _model_coordinate_snapshot(mol_id) if snapshot is None else snapshot
  context_residue_specs = _tq_residue_specs_near_residue_centres(
    mol_id,
    score_residue_specs,
    radius=context_radius,
    max_residues=None,
    snapshot=snapshot,
  )
  selected = set(context_residue_specs)
  score_points = [atom[2] for spec in score_residue_specs for atom in snapshot.get(tuple(spec), {}).get("atoms", [])]
  if score_points and influence_radius > 0:
    centre = tuple(sum(point[axis] for point in score_points) / len(score_points) for axis in range(3))
    bound_sq = (max(math.sqrt(_distance_sq(point, centre)) for point in score_points) + influence_radius)**2
    selected.update(spec for spec, record in snapshot.items()
                    if any(_distance_sq(atom[2], centre) <= bound_sq for atom in record["atoms"]))
  atom_records = []
  for residue_spec in sorted(selected):
    atom_records.extend(_tq_residue_atom_records(mol_id, residue_spec, heavy_atoms_only=True, snapshot=snapshot))
  return atom_records


def _tq_prune_context_atoms(score_atom_records, context_atom_records, sigma, sample_radius):
  influence_radius = float(sample_radius) + float(sigma) * THRESHOLDED_Q_CONTEXT_ATOM_INFLUENCE_SIGMA
  index = _spatial_index(enumerate(context_atom_records), lambda item: item[1]["xyz"], influence_radius)
  keep = {i for atom in score_atom_records for i, _record in _spatial_query(index, atom["xyz"], influence_radius)}
  pruned_context_atoms = [record for i, record in enumerate(context_atom_records) if i in keep]
  return pruned_context_atoms or list(context_atom_records)


def _tq_candidate_grid_indices(score_atom_records, sample_radius, grid_spacing):
  sample_radius_sq = sample_radius * sample_radius
  candidate_indices = set()
  for atom_record in score_atom_records:
    atom_x, atom_y, atom_z = atom_record["xyz"]
    x_index_range = range(
      int(math.floor((atom_x - sample_radius) / grid_spacing)),
      int(math.ceil((atom_x + sample_radius) / grid_spacing)) + 1,
    )
    y_index_range = range(
      int(math.floor((atom_y - sample_radius) / grid_spacing)),
      int(math.ceil((atom_y + sample_radius) / grid_spacing)) + 1,
    )
    z_index_range = range(
      int(math.floor((atom_z - sample_radius) / grid_spacing)),
      int(math.ceil((atom_z + sample_radius) / grid_spacing)) + 1,
    )
    for ix in x_index_range:
      x = ix * grid_spacing
      dx_sq = (x - atom_x) * (x - atom_x)
      if dx_sq > sample_radius_sq:
        continue
      for iy in y_index_range:
        y = iy * grid_spacing
        dxy_sq = dx_sq + (y - atom_y) * (y - atom_y)
        if dxy_sq > sample_radius_sq:
          continue
        for iz in z_index_range:
          z = iz * grid_spacing
          distance_sq = dxy_sq + (z - atom_z) * (z - atom_z)
          if distance_sq > sample_radius_sq:
            continue
          candidate_indices.add((ix, iy, iz))
  return sorted(candidate_indices)


def _tq_numpy_module():
  module = globals().get("_THRESHOLDED_Q_NUMPY", False)
  if module is False:
    try:
      import numpy as module
    except ImportError:
      module = None
    globals()["_THRESHOLDED_Q_NUMPY"] = module
  return module


def _tq_sample_density(map_id, x, y, z):
  try:
    value = float(density_at_point(map_id, x, y, z))
    if not math.isfinite(value):
      raise ValueError("non-finite density")
    return value
  except Exception as error:
    raise RuntimeError(
      f"Map #{map_id} density sampling failed at ({x:.3f}, {y:.3f}, {z:.3f}): {error}"
    ) from error


def _tq_residue_point_records_numpy(np, map_id, score_atoms, context_atoms, sigma, radius, spacing):
  indices = _tq_candidate_grid_indices(score_atoms, radius, spacing)
  coordinates = np.asarray([atom["xyz"] for atom in context_atoms], dtype=float)
  occupancies = np.asarray([atom.get("occupancy", 1.0) for atom in context_atoms], dtype=float)
  score_keys = {atom["atom_key"] for atom in score_atoms}
  owners = np.asarray([atom["atom_key"] in score_keys for atom in context_atoms], dtype=bool)
  # Bound temporary point-by-atom matrices, even in very crowded regions.
  batch_size = max(1, min(512, 262144 // len(context_atoms)))
  records = []
  for start in range(0, len(indices), batch_size):
    points = np.asarray(indices[start:start+batch_size], dtype=float) * spacing
    distances = np.zeros((len(points), len(context_atoms)), dtype=float)
    for axis in range(3):
      delta = points[:, axis, None] - coordinates[None, :, axis]
      distances += delta * delta
    keep = owners[np.argmin(distances, axis=1)]
    points = points[keep]
    distances = distances[keep]
    distances *= -1.0 / (2.0 * sigma * sigma)
    np.exp(distances, out=distances)
    distances *= occupancies
    ideals = distances.sum(axis=1)
    for point, ideal in zip(points, ideals):
      x, y, z = (float(value) for value in point)
      records.append({"experimental_density": _tq_sample_density(map_id, x, y, z),
                      "ideal_density": float(ideal)})
  return records


def _tq_residue_point_records(map_id, score_atom_records, context_atom_records, sigma, sample_radius, grid_spacing):
  if not score_atom_records or not context_atom_records:
    return []
  np = _tq_numpy_module()
  if np is not None:
    return _tq_residue_point_records_numpy(np, map_id, score_atom_records, context_atom_records,
                                            sigma, sample_radius, grid_spacing)

  score_atom_keys = {atom_record["atom_key"] for atom_record in score_atom_records}
  context_atom_data = [
    (*context_atom["xyz"], context_atom["atom_key"])
    for context_atom in context_atom_records
  ]
  occupancies = [float(atom.get("occupancy", 1.0)) for atom in context_atom_records]
  distances = [0.0] * len(context_atom_data)
  gaussian_denominator = 2.0 * sigma * sigma
  exp = math.exp
  point_records = []
  candidate_indices = _tq_candidate_grid_indices(score_atom_records, sample_radius, grid_spacing)

  for ix, iy, iz in candidate_indices:
    x = ix * grid_spacing
    y = iy * grid_spacing
    z = iz * grid_spacing
    # Every candidate is already inside a scoring atom's sampling sphere.
    nearest_context_atom_key = None
    nearest_context_distance_sq = None
    for index, (atom_x, atom_y, atom_z, context_atom_key) in enumerate(context_atom_data):
      dx, dy, dz = x-atom_x, y-atom_y, z-atom_z
      distance_sq = dx*dx + dy*dy + dz*dz
      distances[index] = distance_sq
      if nearest_context_distance_sq is None or distance_sq < nearest_context_distance_sq:
        nearest_context_distance_sq = distance_sq
        nearest_context_atom_key = context_atom_key

    if nearest_context_atom_key is None or nearest_context_atom_key not in score_atom_keys:
      continue

    ideal_density = 0.0
    for distance_sq, occupancy in zip(distances, occupancies):
      ideal_density += occupancy * exp(-distance_sq / gaussian_denominator)

    experimental_density = _tq_sample_density(map_id, x, y, z)

    point_records.append(
      {
        "experimental_density": experimental_density,
        "ideal_density": ideal_density,
      }
    )

  return point_records


def _tq_ideal_density_threshold(point_records, ideal_enclosed_mass_percent):
  positive_densities = sorted(
    (point_record["ideal_density"] for point_record in point_records if point_record["ideal_density"] > 0.0),
    reverse=True,
  )
  if not positive_densities:
    return 0.0
  total_mass = sum(positive_densities)
  target_mass = total_mass * max(0.0, min(100.0, float(ideal_enclosed_mass_percent))) / 100.0
  if target_mass <= 0.0:
    return positive_densities[0] + 1.0e-6
  cumulative_mass = 0.0
  threshold = positive_densities[-1]
  for density_value in positive_densities:
    cumulative_mass += density_value
    threshold = density_value
    if cumulative_mass >= target_mass:
      break
  return threshold


def _tq_residue_region_score(point_records, contour_level, ideal_density_threshold):
  ideal_points = 0
  experimental_points = 0
  overlap_points = 0
  for point_record in point_records:
    ideal_on = point_record["ideal_density"] >= ideal_density_threshold and point_record["ideal_density"] > 0.0
    experimental_on = point_record["experimental_density"] >= contour_level
    if ideal_on:
      ideal_points += 1
    if experimental_on:
      experimental_points += 1
    if ideal_on and experimental_on:
      overlap_points += 1

  denominator = ideal_points + experimental_points
  score = 0.0 if denominator <= 0 else (2.0 * overlap_points) / float(denominator)
  return {
    "score": score,
    "ideal_points": ideal_points,
    "experimental_points": experimental_points,
    "overlap_points": overlap_points,
    "n_points": len(point_records),
  }


def _tq_residue_region_sweep_data(point_records, ideal_density_threshold):
  ideal_flags = [
    (point_record["ideal_density"] >= ideal_density_threshold and point_record["ideal_density"] > 0.0)
    for point_record in point_records
  ]
  ideal_points = sum(1 for ideal_on in ideal_flags if ideal_on)
  sorted_pairs = sorted(
    zip((point_record["experimental_density"] for point_record in point_records), ideal_flags),
    key=lambda pair: pair[0],
  )
  experimental_densities = [pair[0] for pair in sorted_pairs]
  overlap_suffix = [0] * (len(sorted_pairs) + 1)
  for index in range(len(sorted_pairs) - 1, -1, -1):
    overlap_suffix[index] = overlap_suffix[index + 1] + (1 if sorted_pairs[index][1] else 0)
  return {
    "ideal_points": ideal_points,
    "n_points": len(point_records),
    "experimental_densities": experimental_densities,
    "overlap_suffix": overlap_suffix,
  }


def _tq_residue_region_score_from_sweep_data(sweep_data, contour_level):
  experimental_densities = sweep_data["experimental_densities"]
  idx = bisect_left(experimental_densities, contour_level)
  experimental_points = len(experimental_densities) - idx
  overlap_points = sweep_data["overlap_suffix"][idx]
  ideal_points = sweep_data["ideal_points"]
  denominator = ideal_points + experimental_points
  score = 0.0 if denominator <= 0 else (2.0 * overlap_points) / float(denominator)
  return {
    "score": score,
    "ideal_points": ideal_points,
    "experimental_points": experimental_points,
    "overlap_points": overlap_points,
    "n_points": sweep_data["n_points"],
  }


def _tq_prepare_scored_residues(
    map_id,
    mol_id,
    residue_specs,
    sigma=THRESHOLDED_Q_GAUSSIAN_SIGMA,
    ideal_enclosed_mass_percent=THRESHOLDED_Q_IDEAL_ENCLOSED_MASS_PERCENT,
    grid_spacing=THRESHOLDED_Q_GRID_SPACING,
    sample_shell_padding=THRESHOLDED_Q_SAMPLE_SHELL_PADDING,
    context_radius=None,
    backbone_only=False,
    snapshot=None,
):
  try:
    sigma = float(sigma)
    grid_spacing = float(grid_spacing)
    sample_shell_padding = float(sample_shell_padding)
    ideal_enclosed_mass_percent = float(ideal_enclosed_mass_percent)
  except Exception:
    return None
  if not all(math.isfinite(value) for value in (sigma, grid_spacing, sample_shell_padding, ideal_enclosed_mass_percent)) or sigma <= 0.0 or grid_spacing <= 0.0 or sample_shell_padding < 0.0 or not 0 <= ideal_enclosed_mass_percent <= 100:
    return None

  score_residue_specs = []
  for residue_spec in residue_specs:
    parsed = _coot_residue_spec_from_spec(residue_spec)
    if parsed is not None and parsed not in score_residue_specs:
      score_residue_specs.append(parsed)
  if not score_residue_specs:
    return None

  support_radius = _tq_support_radius(sigma, float(ideal_enclosed_mass_percent))
  sample_radius = sigma * THRESHOLDED_Q_SAMPLE_RADIUS_SIGMA + sample_shell_padding
  if context_radius is None:
    context_radius = THRESHOLDED_Q_LOCAL_RADIUS
  context_radius = max(float(context_radius), float(sample_radius))
  snapshot = _model_coordinate_snapshot(mol_id) if snapshot is None else snapshot
  influence_radius = sample_radius + sigma * THRESHOLDED_Q_CONTEXT_ATOM_INFLUENCE_SIGMA
  context_atom_records = _tq_context_atom_records(mol_id, score_residue_specs, context_radius,
                                                snapshot=snapshot, influence_radius=influence_radius)
  if not context_atom_records:
    return None

  prepared_residues = []
  for residue_spec in score_residue_specs:
    score_atom_records = _tq_score_atom_records(
      _tq_residue_atom_records(mol_id, residue_spec, heavy_atoms_only=True, snapshot=snapshot),
      backbone_only=backbone_only,
    )
    if not score_atom_records:
      continue
    residue_context_atoms = _tq_prune_context_atoms(score_atom_records, context_atom_records, sigma, sample_radius)
    point_records = _tq_residue_point_records(
      map_id,
      score_atom_records,
      residue_context_atoms,
      sigma,
      sample_radius,
      grid_spacing,
    )
    if not point_records:
      continue
    residue_name_here = score_atom_records[0]["residue_name"]
    ideal_density_threshold = _tq_ideal_density_threshold(point_records, ideal_enclosed_mass_percent)
    prepared_residues.append(
      {
        "residue_spec": residue_spec,
        "residue_name": residue_name_here,
        "label": _tq_residue_label(residue_spec[0], residue_spec[1], residue_spec[2], residue_name_here),
        "n_atoms": len(score_atom_records),
        "n_points": len(point_records),
        "ideal_density_threshold": ideal_density_threshold,
        "region_sweep_data": _tq_residue_region_sweep_data(point_records, ideal_density_threshold),
      }
    )

  if not prepared_residues:
    return None

  return {
    "map_id": map_id,
    "mol_id": mol_id,
    "score_residue_specs": score_residue_specs,
    "support_radius": support_radius,
    "sample_radius": sample_radius,
    "sigma": float(sigma),
    "ideal_enclosed_mass_percent": float(ideal_enclosed_mass_percent),
    "prepared_residues": prepared_residues,
  }


def _tq_prepare_scored_residues_or_status(*args, quiet=False, **kwargs):
  try:
    prepared = _tq_prepare_scored_residues(*args, **kwargs)
    if prepared is None and not quiet:
      _tq_status("Thresholded Q-residue: no scorable atoms or invalid sampling parameters")
    return prepared
  except (RuntimeError, ValueError, TypeError, OverflowError) as error:
    # Never optimize against a partially sampled map or fabricated zero density.
    if not quiet:
      _tq_status(f"Thresholded Q-residue: calculation aborted. {error}")
    return None


def thresholded_q_residues(
    map_id,
    mol_id,
    residue_specs,
    contour_level=None,
    sigma=THRESHOLDED_Q_GAUSSIAN_SIGMA,
    ideal_enclosed_mass_percent=THRESHOLDED_Q_IDEAL_ENCLOSED_MASS_PERCENT,
    grid_spacing=THRESHOLDED_Q_GRID_SPACING,
    sample_shell_padding=THRESHOLDED_Q_SAMPLE_SHELL_PADDING,
    context_radius=None,
    backbone_only=False,
    quiet=False,
):
  map_id = _tq_map_id_or_default(map_id)
  try:
    sigma = float(sigma)
    grid_spacing = float(grid_spacing)
    sample_shell_padding = float(sample_shell_padding)
    if context_radius is not None:
      context_radius = float(context_radius)
  except Exception:
    message = "Thresholded Q-residue: invalid numeric parameter"
    if not quiet:
      _tq_status(message)
    return None
  if sigma <= 0.0 or grid_spacing <= 0.0 or sample_shell_padding < 0.0:
    message = "Thresholded Q-residue: sigma/grid spacing must be > 0 and shell padding must be >= 0"
    if not quiet:
      _tq_status(message)
    return None
  if not _tq_valid_map_molecule(map_id):
    message = f"Thresholded Q-residue: invalid map #{map_id}"
    if not quiet:
      _tq_status(message)
    return None
  if not _tq_valid_model_molecule(mol_id):
    message = f"Thresholded Q-residue: invalid model molecule #{mol_id}"
    if not quiet:
      _tq_status(message)
    return None
  if contour_level is None:
    try:
      contour_level = float(get_contour_level_absolute(map_id))
    except Exception:
      contour_level = None
  if contour_level is None:
    contour_level = 0.0
  try:
    contour_level = float(contour_level)
    if not math.isfinite(contour_level):
      raise ValueError("non-finite contour")
  except (TypeError, ValueError, OverflowError):
    if not quiet:
      _tq_status("Thresholded Q-residue: invalid contour level")
    return None

  prepared = _tq_prepare_scored_residues_or_status(
    map_id,
    mol_id,
    residue_specs,
    sigma=sigma,
    ideal_enclosed_mass_percent=ideal_enclosed_mass_percent,
    grid_spacing=grid_spacing,
    sample_shell_padding=sample_shell_padding,
    context_radius=context_radius,
    backbone_only=backbone_only,
    quiet=quiet,
  )
  if prepared is None:
    return None

  residue_scores = []
  for prepared_residue in prepared["prepared_residues"]:
    region_score = _tq_residue_region_score_from_sweep_data(
      prepared_residue["region_sweep_data"],
      contour_level,
    )
    residue_scores.append(
      {
        "residue_spec": prepared_residue["residue_spec"],
        "residue_name": prepared_residue["residue_name"],
        "label": prepared_residue["label"],
        "score": region_score["score"],
        "n_atoms": prepared_residue["n_atoms"],
        "n_points": region_score["n_points"],
        "ideal_points": region_score["ideal_points"],
        "experimental_points": region_score["experimental_points"],
        "overlap_points": region_score["overlap_points"],
        "ideal_density_threshold": prepared_residue["ideal_density_threshold"],
      }
    )

  overall_score = 0.0
  if residue_scores:
    overall_score = sum(residue_score["score"] for residue_score in residue_scores) / float(len(residue_scores))

  result = {
    "map_id": map_id,
    "mol_id": mol_id,
    "contour_level": contour_level,
    "score": overall_score,
    "residue_scores": residue_scores,
    "support_radius": prepared["support_radius"],
    "sample_radius": prepared["sample_radius"],
    "sigma": prepared["sigma"],
    "ideal_enclosed_mass_percent": prepared["ideal_enclosed_mass_percent"],
  }

  if not quiet:
    labels = ", ".join(residue_score["label"] for residue_score in residue_scores[:4])
    if len(residue_scores) > 4:
      labels += ", ..."
    _tq_status(
      "Thresholded Q-residue at level {0:.4f}: score {1:.3f} for {2} residue(s) in {3} [{4}]".format(
        contour_level,
        overall_score,
        len(residue_scores),
        _molecule_id_name_label(mol_id),
        labels or "no labels",
      )
    )
  return result


def thresholded_q_residue(map_id, mol_id, residue_spec, **kwargs):
  result = thresholded_q_residues(map_id, mol_id, [residue_spec], **kwargs)
  if result is None:
    return None
  if result["residue_scores"]:
    result["residue_score"] = result["residue_scores"][0]
  else:
    result["residue_score"] = None
  return result


def thresholded_q_active_residue(
    map_id=None,
    ideal_enclosed_mass_percent=THRESHOLDED_Q_IDEAL_ENCLOSED_MASS_PERCENT,
    **kwargs,
):
  target = _tq_active_residue_target()
  if target is None:
    message = "Thresholded Q-residue: no active residue"
    _tq_status(message)
    return None
  mol_id, chain_id, resno, ins_code = target
  return thresholded_q_residue(
    _tq_map_id_or_default(map_id),
    mol_id,
    (chain_id, resno, ins_code),
    ideal_enclosed_mass_percent=ideal_enclosed_mass_percent,
    **kwargs,
  )


def optimize_local_threshold_by_thresholded_q(
    map_id,
    mol_id,
    residue_specs,
    contour_min=None,
    contour_max=None,
    contour_min_sigma=None,
    contour_max_sigma=None,
    n_steps=THRESHOLDED_Q_SWEEP_STEPS,
    sigma=THRESHOLDED_Q_GAUSSIAN_SIGMA,
    ideal_enclosed_mass_percent=THRESHOLDED_Q_IDEAL_ENCLOSED_MASS_PERCENT,
    grid_spacing=THRESHOLDED_Q_GRID_SPACING,
    sample_shell_padding=THRESHOLDED_Q_SAMPLE_SHELL_PADDING,
    context_radius=None,
    backbone_only=False,
    apply_contour=False,
    quiet=False,
    snapshot=None,
):
  map_id = _tq_map_id_or_default(map_id)
  try:
    sigma = float(sigma)
    grid_spacing = float(grid_spacing)
    sample_shell_padding = float(sample_shell_padding)
    n_steps = int(n_steps)
    ideal_enclosed_mass_percent = float(ideal_enclosed_mass_percent)
    if not all(math.isfinite(value) for value in (sigma, grid_spacing, sample_shell_padding, ideal_enclosed_mass_percent)):
      raise ValueError("Non-finite sampling parameter")
    if not 0.0 < ideal_enclosed_mass_percent <= 100.0:
      raise ValueError("Enclosed mass must be in (0, 100]")
    contour_min = float(contour_min) if contour_min is not None else None
    contour_max = float(contour_max) if contour_max is not None else None
    contour_min_sigma = float(contour_min_sigma) if contour_min_sigma is not None else None
    contour_max_sigma = float(contour_max_sigma) if contour_max_sigma is not None else None
    if not all(math.isfinite(value) for value in (contour_min, contour_max, contour_min_sigma, contour_max_sigma) if value is not None):
      raise ValueError("Non-finite contour bound")
    if context_radius is not None:
      context_radius = float(context_radius)
      if not math.isfinite(context_radius) or context_radius <= 0:
        raise ValueError("Invalid context radius")
  except Exception:
    message = "Thresholded Q-residue optimizer: invalid numeric parameter"
    if not quiet:
      _tq_status(message)
    return None
  if sigma <= 0.0 or grid_spacing <= 0.0 or sample_shell_padding < 0.0:
    message = "Thresholded Q-residue optimizer: sigma/grid spacing must be > 0 and shell padding must be >= 0"
    if not quiet:
      _tq_status(message)
    return None
  if not _tq_valid_map_molecule(map_id):
    message = f"Thresholded Q-residue optimizer: invalid map #{map_id}"
    if not quiet:
      _tq_status(message)
    return None
  if not _tq_valid_model_molecule(mol_id):
    message = f"Thresholded Q-residue optimizer: invalid model molecule #{mol_id}"
    if not quiet:
      _tq_status(message)
    return None
  try:
    current_contour = float(get_contour_level_absolute(map_id))
  except Exception:
    current_contour = 0.0
  try:
    current_contour_sigma = float(get_contour_level_in_sigma(map_id))
  except Exception:
    current_contour_sigma = None
  map_sigma = _tq_map_sigma(map_id)

  if contour_min is None and contour_max is None and map_sigma is not None:
    if contour_min_sigma is None:
      contour_min_sigma = THRESHOLDED_Q_SWEEP_MIN_SIGMA
    if contour_max_sigma is None:
      contour_max_sigma = THRESHOLDED_Q_SWEEP_MAX_SIGMA
    contour_min = float(contour_min_sigma) * map_sigma
    contour_max = float(contour_max_sigma) * map_sigma
  else:
    if contour_min is None:
      contour_min = float(contour_min_sigma) * map_sigma if contour_min_sigma is not None and map_sigma is not None else current_contour
    if contour_max is None:
      contour_max = float(contour_max_sigma) * map_sigma if contour_max_sigma is not None and map_sigma is not None else current_contour
  contour_min = float(contour_min)
  contour_max = float(contour_max)
  if contour_max < contour_min:
    contour_min, contour_max = contour_max, contour_min
  n_steps = max(2, n_steps)

  prepared = _tq_prepare_scored_residues_or_status(
    map_id,
    mol_id,
    residue_specs,
    sigma=sigma,
    ideal_enclosed_mass_percent=ideal_enclosed_mass_percent,
    grid_spacing=grid_spacing,
    sample_shell_padding=sample_shell_padding,
    context_radius=context_radius,
    backbone_only=backbone_only,
    quiet=quiet,
    snapshot=snapshot,
  )
  if prepared is None:
    return None

  contour_levels = [contour_min] if contour_max == contour_min else [
    contour_min + (contour_max - contour_min) * (index / float(n_steps - 1))
    for index in range(n_steps)
  ]

  sweep_results = []
  best_result = None
  for contour_level in contour_levels:
    residue_scores = []
    for prepared_residue in prepared["prepared_residues"]:
      region_score = _tq_residue_region_score_from_sweep_data(
        prepared_residue["region_sweep_data"],
        contour_level,
      )
      residue_scores.append(region_score["score"])
    mean_score = sum(residue_scores) / float(len(residue_scores)) if residue_scores else 0.0
    sweep_result = {
      "contour_level": contour_level,
      "score": mean_score,
      "n_residues": len(residue_scores),
    }
    if map_sigma is not None and map_sigma > 0.0:
      sweep_result["contour_sigma"] = contour_level / map_sigma
    sweep_results.append(sweep_result)
    if best_result is None:
      best_result = sweep_result
    elif mean_score > best_result["score"] + 1.0e-9:
      best_result = sweep_result
    elif abs(mean_score - best_result["score"]) <= 1.0e-9 and contour_level > best_result["contour_level"]:
      best_result = sweep_result

  result = {
    "map_id": map_id,
    "mol_id": mol_id,
    "current_contour": current_contour,
    "current_contour_sigma": current_contour_sigma,
    "map_sigma": map_sigma,
    "best_contour_level": best_result["contour_level"],
    "best_score": best_result["score"],
    "sweep_results": sweep_results,
    "residue_specs": prepared["score_residue_specs"],
    "support_radius": prepared["support_radius"],
    "sample_radius": prepared["sample_radius"],
    "ideal_enclosed_mass_percent": prepared["ideal_enclosed_mass_percent"],
  }
  if map_sigma is not None and map_sigma > 0.0:
    result["best_contour_sigma"] = result["best_contour_level"] / map_sigma

  if not quiet:
    if "best_contour_sigma" in result:
      _tq_status(
        "Thresholded Q-residue local optimum: level {0:.4f}, {1:.2f} sigma (score {2:.3f}) for {3} residue(s) in {4}".format(
          result["best_contour_level"],
          result["best_contour_sigma"],
          result["best_score"],
          len(result["residue_specs"]),
          _molecule_id_name_label(mol_id),
        )
      )
    else:
      _tq_status(
        "Thresholded Q-residue local optimum: level {0:.4f} (score {1:.3f}) for {2} residue(s) in {3}".format(
          result["best_contour_level"],
          result["best_score"],
          len(result["residue_specs"]),
          _molecule_id_name_label(mol_id),
        )
      )

  if apply_contour:
    try:
      set_contour_level_absolute(map_id, result["best_contour_level"])
      _tq_post_map_level_status(map_id, "set to")
    except Exception as error:
      _tq_status(f"Thresholded Q-residue optimizer: failed to apply contour ({error!r})")
  return result


def optimize_local_threshold_by_thresholded_q_current_view(
    map_id=None,
    mol_id=None,
    residue_specs=None,
    radius=THRESHOLDED_Q_LOCAL_RADIUS,
    max_residues=THRESHOLDED_Q_LOCAL_MAX_RESIDUES,
    apply_contour=True,
    **kwargs,
):
  map_id = _tq_map_id_or_default(map_id)
  try:
    radius = float(radius)
    max_residues = int(max_residues)
  except Exception:
    _tq_status("Thresholded Q-residue optimizer: invalid radius/max_residues")
    return None
  if not math.isfinite(radius) or radius <= 0.0:
    _tq_status("Thresholded Q-residue optimizer: radius must be > 0")
    return None
  if max_residues <= 0:
    _tq_status("Thresholded Q-residue optimizer: max_residues must be > 0")
    return None
  if residue_specs is None:
    target = _tq_active_residue_target()
    if target is None:
      _tq_status("Thresholded Q-residue optimizer: no active residue")
      return None
    mol_id = target[0]
    residue_specs = [(target[1], target[2], target[3])]
  elif mol_id is None:
    target = _tq_active_residue_target()
    if target is not None:
      mol_id = target[0]
  if mol_id is None:
    _tq_status("Thresholded Q-residue optimizer: no target molecule")
    return None

  try:
    snapshot = _model_coordinate_snapshot(mol_id)
    score_residue_specs = _tq_residue_specs_near_residue_centres(
      mol_id, residue_specs, radius=radius, max_residues=max_residues, snapshot=snapshot)
  except (RuntimeError, ValueError, TypeError) as error:
    _tq_status(f"Thresholded Q-residue optimizer: calculation aborted. {error}")
    return None
  if not score_residue_specs:
    _tq_status("Thresholded Q-residue optimizer: no local residues to score")
    return None
  context_radius = kwargs.pop("context_radius", radius + THRESHOLDED_Q_CONTEXT_RADIUS_PADDING)
  if context_radius is None:
    context_radius = radius + THRESHOLDED_Q_CONTEXT_RADIUS_PADDING
  try:
    context_radius = float(context_radius)
    if not math.isfinite(context_radius) or context_radius <= 0:
      raise ValueError("context radius must be finite and positive")
  except (ValueError, TypeError) as error:
    _tq_status(f"Thresholded Q-residue optimizer: {error}")
    return None
  return optimize_local_threshold_by_thresholded_q(
    map_id,
    mol_id,
    score_residue_specs,
    context_radius=max(radius, context_radius),
    apply_contour=apply_contour,
    snapshot=snapshot,
    **kwargs,
  )


def _step_map_sigma(mol_id, delta):
  """Adjust contour sigma in coarse 0.1 steps within the usual trimmings range."""
  current_level = get_contour_level_in_sigma(mol_id)
  if 0.1 <= current_level <= 10.0:
    new_level = current_level + delta
  elif current_level < 0.1:
    new_level = 0.1
  else:
    new_level = 10.0
  set_contour_level_in_sigma(mol_id, new_level)
  _post_map_level_status(mol_id, "set to")


def step_map_coarse_up(mol_id):
  _step_map_sigma(mol_id, 0.1)


def step_map_coarse_down(mol_id):
  _step_map_sigma(mol_id, -0.1)


def set_current_map_sigma(level):
  map_id = _scrollable_map_or_status()
  if map_id is None:
    return None
  set_contour_level_in_sigma(map_id, level)
  _post_map_level_status(map_id, "set to")


def step_current_map_coarse_up():
  map_id = _scrollable_map_or_status()
  if map_id is None:
    return None
  step_map_coarse_up(map_id)


def step_current_map_coarse_down():
  map_id = _scrollable_map_or_status()
  if map_id is None:
    return None
  step_map_coarse_down(map_id)


def _map_brightness_scale_would_clip(map_id, scale_factor):
  """Return True if a brighten/darken step would trigger Coot's colour clamp."""
  try:
    current_colour = map_colour_components_py(map_id)
  except Exception:
    return False
  if not isinstance(current_colour, (list, tuple)):
    return False
  try:
    for component in current_colour:
      proposed = float(component) * float(scale_factor)
      if proposed < MAP_BRIGHTNESS_MIN_COMPONENT or proposed > MAP_BRIGHTNESS_MAX_COMPONENT:
        return True
  except Exception:
    return False
  return False


def _scale_scrollable_map_brightness(scale_factor, status_message, clipped_message):
  """Brighten or darken only the current scrollable map."""
  map_id = _scrollable_map_or_status()
  if map_id is None:
    return None
  if _map_brightness_scale_would_clip(map_id, scale_factor):
    add_status_bar_text(clipped_message)
    return None
  coot_utils.brighten_map(map_id, scale_factor)
  add_status_bar_text(status_message)
  return map_id


def brighten_scrollable_map():
  return _scale_scrollable_map_brightness(
    MAP_BRIGHTEN_SCALE_FACTOR,
    "Brightened current map",
    "Skipped brightening current map to avoid hitting map-colour limits",
  )


def darken_scrollable_map():
  return _scale_scrollable_map_brightness(
    MAP_DARKEN_SCALE_FACTOR,
    "Darkened current map",
    "Skipped darkening current map to avoid hitting map-colour limits",
  )


def start_measure_distance():
  add_status_bar_text("Pick 2 atoms to measure distance")
  do_distance_define()


def _map_looks_em_like(map_id):
  if not is_valid_map_molecule(map_id):
    return False
  try:
    cell = map_cell(map_id)
  except Exception:
    cell = None
  if cell and len(cell) >= 6:
    alpha, beta, gamma = cell[3], cell[4], cell[5]
    if (
      abs(alpha - 90.0) < 0.5
      and abs(beta - 90.0) < 0.5
      and abs(gamma - 90.0) < 0.5
    ):
      return True
  for path in _candidate_map_file_paths(map_id):
    if path.lower().endswith((".mrc", ".map", ".ccp4")):
      return True
  return False


def _candidate_map_file_paths(map_id):
  name = molecule_name(map_id)
  if not name:
    return []
  candidates = [name]
  if not os.path.isabs(name):
    candidates.append(os.path.join(os.getcwd(), name))
  return [path for path in candidates if os.path.isfile(path)]


def _parse_ccp4_mrc_header(path):
  with open(path, "rb") as handle:
    header = handle.read(1024)
  if len(header) < 76:
    return None
  for endian in ("<", ">"):
    ints_10 = struct.unpack(endian + "10i", header[:40])
    mx, my, mz = ints_10[7], ints_10[8], ints_10[9]
    floats_6 = struct.unpack(endian + "6f", header[40:64])
    a, b, c, alpha, beta, gamma = floats_6
    axes = struct.unpack(endian + "3i", header[64:76])
    if (
      mx > 0 and my > 0 and mz > 0
      and a > 0.0 and b > 0.0 and c > 0.0
      and 0.0 < alpha <= 180.0
      and 0.0 < beta <= 180.0
      and 0.0 < gamma <= 180.0
      and all(axis in (1, 2, 3) for axis in axes)
    ):
      return {
        "mx": mx,
        "my": my,
        "mz": mz,
        "a": a,
        "b": b,
        "c": c,
      }
  return None


def _resample_factor_for_target_pixel_size(map_id, target_pixel_size):
  for path in _candidate_map_file_paths(map_id):
    header = _parse_ccp4_mrc_header(path)
    if not header:
      continue
    sampling_a = header["a"] / float(header["mx"])
    sampling_b = header["b"] / float(header["my"])
    sampling_c = header["c"] / float(header["mz"])
    current_pixel_size = max(sampling_a, sampling_b, sampling_c)
    return current_pixel_size / target_pixel_size, current_pixel_size
  return None, None


def _resample_plan_for_target_pixel_size(map_id, target_pixel_size):
  """Estimate the output grid that resampling would create for a target pixel size."""
  for path in _candidate_map_file_paths(map_id):
    header = _parse_ccp4_mrc_header(path)
    if not header:
      continue
    sampling_a = header["a"] / float(header["mx"])
    sampling_b = header["b"] / float(header["my"])
    sampling_c = header["c"] / float(header["mz"])
    current_pixel_size = max(sampling_a, sampling_b, sampling_c)
    resample_factor = current_pixel_size / target_pixel_size
    return {
      "resample_factor": resample_factor,
      "current_pixel_size": current_pixel_size,
      "new_grid": (
        int(math.ceil(header["mx"] * resample_factor)),
        int(math.ceil(header["my"] * resample_factor)),
        int(math.ceil(header["mz"] * resample_factor)),
      ),
    }
  return None


def style_resampled_em_map(map_id):
  """Apply the standard ambient EM-map appearance used by the restyle shortcut."""
  set_draw_map_standard_lines(map_id, 1)
  set_draw_solid_density_surface(map_id, 0)
  if not map_is_difference_map(map_id):
    set_map_colour(map_id, *EM_REFINED_MAP_COLOUR)
  for name, values in (("set_map_material_ambient", EM_MAP_AMBIENT),
                       ("set_map_material_diffuse", EM_MAP_DIFFUSE)):
    setter = getattr(coot, name, None)
    if callable(setter):
      setter(map_id, *values)
  set_map_material_specular(map_id, *EM_MAP_SPECULAR)
  set_map_fresnel_settings(map_id, *EM_MAP_FRESNEL)
  set_solid_density_surface_opacity(map_id, 1.0)
  MAP_LOCAL_APPEARANCE_STATE[map_id] = {
    "map_colour": list(map_colour_components_py(map_id)),
    "ambient": EM_MAP_AMBIENT,
    "diffuse": EM_MAP_DIFFUSE,
    "specular_strength": EM_MAP_SPECULAR[0],
    "shininess": EM_MAP_SPECULAR[1],
    "fresnel": EM_MAP_FRESNEL,
    "opacity": 1.0,
  }


def _map_appearance_map_label():
  label = Gtk.Label(label="")
  label.set_xalign(0.0)
  label.set_hexpand(True)
  label.set_ellipsize(_coot_gui_repository_module("Pango").EllipsizeMode.MIDDLE)
  label.set_max_width_chars(44)
  return label


def _update_map_appearance_map_label(label, map_id):
  if map_id is None:
    label.set_text("No active map")
    label.set_tooltip_text(None)
    return
  try:
    name = molecule_name(map_id)
  except Exception:
    name = "map"
  text = f"Active map: {name} ({map_id})"
  label.set_text(text)
  label.set_tooltip_text(text)


def show_active_map_brightness_slider():
  """Open a non-modal slider that live-updates the current active map's ambient brightness."""
  initial_map = _scrollable_map_or_status()
  if initial_map is None:
    return None

  window = Gtk.Window()
  window.set_title("Map brightness")
  window.set_default_size(420, -1)

  vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=10)
  vbox.set_margin_top(12)
  vbox.set_margin_bottom(12)
  vbox.set_margin_start(12)
  vbox.set_margin_end(12)

  map_label = _map_appearance_map_label()
  value_label = Gtk.Label(label="100%")
  value_label.set_width_chars(5)

  slider_row = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=10)
  adjustment = Gtk.Adjustment(value=100.0, lower=0.0, upper=200.0, step_increment=1.0, page_increment=10.0)
  slider = Gtk.Scale(orientation=Gtk.Orientation.HORIZONTAL, adjustment=adjustment)
  slider.set_hexpand(True)
  slider.set_draw_value(False)
  slider.set_digits(0)

  button_row = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
  button_row.set_halign(Gtk.Align.END)
  reset_button = Gtk.Button(label="Reset to 100%")
  close_button = Gtk.Button(label="Close")

  def current_map():
    # Resolve this afresh for every update so changing the scroll-wheel/refinement
    # map while the window is open redirects the control to that map.
    map_id = scroll_wheel_map()
    if map_id != -1 and map_id in map_molecule_list():
      return map_id
    map_id = imol_refinement_map()
    if map_id != -1 and map_id in map_molecule_list():
      return map_id
    displayed = [imol for imol in map_molecule_list() if map_is_displayed(imol)]
    return displayed[0] if displayed else None

  def update_map_label(map_id):
    _update_map_appearance_map_label(map_label, map_id)

  def apply_brightness(scale):
    percent = float(scale.get_value())
    value_label.set_text("{:.0f}%".format(percent))
    map_id = current_map()
    update_map_label(map_id)
    if map_id is None:
      return
    brightness = percent / 100.0
    try:
      set_map_material_ambient(map_id, brightness, brightness, brightness, 1.0)
      state = MAP_LOCAL_APPEARANCE_STATE.setdefault(map_id, {})
      state["ambient"] = (brightness, brightness, brightness, 1.0)
    except Exception as error:
      print("WARNING:: could not set map brightness: {}".format(error))

  def reset_brightness(button):
    slider.set_value(100.0)

  slider.connect("value-changed", apply_brightness)
  reset_button.connect("clicked", reset_brightness)
  close_button.connect("clicked", lambda button: window.close())

  slider_row.append(slider)
  slider_row.append(value_label)
  button_row.append(reset_button)
  button_row.append(close_button)
  vbox.append(map_label)
  vbox.append(slider_row)
  vbox.append(button_row)
  window.set_child(vbox)
  update_map_label(initial_map)
  window.present()
  return window


def show_active_map_colour_slider():
  """Open a non-modal hue slider that live-updates the current active map colour."""
  import colorsys

  initial_map = _scrollable_map_or_status()
  if initial_map is None:
    return None

  window = Gtk.Window()
  window.set_title("Map colour")
  window.set_default_size(420, -1)

  vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=10)
  vbox.set_margin_top(12)
  vbox.set_margin_bottom(12)
  vbox.set_margin_start(12)
  vbox.set_margin_end(12)

  map_label = _map_appearance_map_label()
  value_label = Gtk.Label(label="")
  value_label.set_width_chars(5)

  slider_row = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=10)
  adjustment = Gtk.Adjustment(value=0.0, lower=0.0, upper=360.0, step_increment=1.0, page_increment=15.0)
  slider = Gtk.Scale(orientation=Gtk.Orientation.HORIZONTAL, adjustment=adjustment)
  slider.set_hexpand(True)
  slider.set_draw_value(False)
  slider.set_digits(0)

  button_row = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
  button_row.set_halign(Gtk.Align.END)
  reset_button = Gtk.Button(label="Reset colour")
  close_button = Gtk.Button(label="Close")

  state = {"map_id": None, "saturation": 1.0, "value": 1.0, "reset_colour": None, "updating": False}

  def current_map():
    map_id = scroll_wheel_map()
    if map_id != -1 and map_id in map_molecule_list():
      return map_id
    map_id = imol_refinement_map()
    if map_id != -1 and map_id in map_molecule_list():
      return map_id
    displayed = [imol for imol in map_molecule_list() if map_is_displayed(imol)]
    return displayed[0] if displayed else None

  def update_map_label(map_id):
    _update_map_appearance_map_label(map_label, map_id)

  def sync_to_map(map_id):
    if map_id is None or map_id == state["map_id"]:
      return
    try:
      colour = tuple(float(x) for x in map_colour_components_py(map_id)[:3])
      hue, saturation, value = colorsys.rgb_to_hsv(*colour)
    except Exception as error:
      print("WARNING:: could not read map colour: {}".format(error))
      return
    state["map_id"] = map_id
    state["saturation"] = saturation
    state["value"] = value
    state["reset_colour"] = colour
    state["updating"] = True
    slider.set_value(hue * 360.0)
    state["updating"] = False
    value_label.set_text("{:.0f} deg".format(hue * 360.0))

  def apply_colour(scale):
    if state["updating"]:
      return
    map_id = current_map()
    update_map_label(map_id)
    if map_id is None:
      return
    sync_to_map(map_id)
    degrees = float(scale.get_value())
    value_label.set_text("{:.0f} deg".format(degrees))
    red, green, blue = colorsys.hsv_to_rgb(degrees / 360.0, state["saturation"], state["value"])
    try:
      set_map_colour(map_id, red, green, blue)
      local_state = MAP_LOCAL_APPEARANCE_STATE.setdefault(map_id, {})
      local_state["map_colour"] = [red, green, blue]
    except Exception as error:
      print("WARNING:: could not set map colour: {}".format(error))

  def reset_colour(button):
    map_id = current_map()
    update_map_label(map_id)
    if map_id is None:
      return
    sync_to_map(map_id)
    colour = state["reset_colour"]
    if colour is None:
      return
    try:
      set_map_colour(map_id, *colour)
      hue, saturation, value = colorsys.rgb_to_hsv(*colour)
      state["saturation"] = saturation
      state["value"] = value
      state["updating"] = True
      slider.set_value(hue * 360.0)
      state["updating"] = False
      value_label.set_text("{:.0f} deg".format(hue * 360.0))
      local_state = MAP_LOCAL_APPEARANCE_STATE.setdefault(map_id, {})
      local_state["map_colour"] = list(colour)
    except Exception as error:
      print("WARNING:: could not reset map colour: {}".format(error))

  slider.connect("value-changed", apply_colour)
  reset_button.connect("clicked", reset_colour)
  close_button.connect("clicked", lambda button: window.close())

  slider_row.append(slider)
  slider_row.append(value_label)
  button_row.append(reset_button)
  button_row.append(close_button)
  vbox.append(map_label)
  vbox.append(slider_row)
  vbox.append(button_row)
  window.set_child(vbox)
  update_map_label(initial_map)
  sync_to_map(initial_map)
  window.present()
  return window


def _restyle_active_map_without_resampling(map_id, contour_level_sigma, status_message):
  style_resampled_em_map(map_id)
  if isinstance(contour_level_sigma, (int, float)):
    set_contour_level_in_sigma(map_id, contour_level_sigma)
  set_map_displayed(map_id, 1)
  set_scroll_wheel_map(map_id)
  set_scrollable_map(map_id)
  add_status_bar_text(status_message)
  return map_id


def resample_active_map_for_em_half_angstrom(force=False, allow_large_output=False):
  map_id = _scrollable_map_or_status()
  if map_id is None:
    return None
  if map_is_difference_map(map_id):
    add_status_bar_text("Scroll-wheel map is a difference map; resample/restyle skipped")
    return None
  original_refinement_map = imol_refinement_map()
  old_map_was_refinement_map = (original_refinement_map == map_id)
  contour_level_sigma = get_contour_level_in_sigma(map_id)
  if not force and not _map_looks_em_like(map_id):
    return _restyle_active_map_without_resampling(
      map_id,
      contour_level_sigma,
      "Map is not an EM map",
    )
  resample_plan = _resample_plan_for_target_pixel_size(map_id, EM_TARGET_PIXEL_SIZE)
  if resample_plan is None:
    return _restyle_active_map_without_resampling(
      map_id,
      contour_level_sigma,
      "Grid spacing could not be determined",
    )
  resample_factor = resample_plan["resample_factor"]
  current_pixel_size = resample_plan["current_pixel_size"]
  if current_pixel_size <= EM_TARGET_PIXEL_SIZE:
    style_resampled_em_map(map_id)
    set_contour_level_in_sigma(map_id, contour_level_sigma)
    set_map_displayed(map_id, 1)
    set_scroll_wheel_map(map_id)
    set_scrollable_map(map_id)
    if old_map_was_refinement_map:
      try:
        set_imol_refinement_map(map_id)
      except Exception:
        pass
    add_status_bar_text("Map already finer than 0.5 A/pixel; restyled without resampling")
    return map_id
  new_grid = resample_plan["new_grid"]
  max_grid_dimension = max(new_grid)
  if (
    not allow_large_output
    and max_grid_dimension > EM_RESAMPLE_CONFIRM_MAX_GRID_DIMENSION
  ):
    generic_confirm_dialog(
      "Resample/restyle current map",
      "Resampling would generate a large map, which may take a while.\n\n"
      "Estimated grid: {0} x {1} x {2} (max dimension {3} > {4}).\n\n"
      "What do you want to do?".format(
        new_grid[0],
        new_grid[1],
        new_grid[2],
        max_grid_dimension,
        EM_RESAMPLE_CONFIRM_MAX_GRID_DIMENSION,
      ),
      "Just restyle - skip resampling",
      lambda: _restyle_active_map_without_resampling(
        map_id,
        contour_level_sigma,
        "Skipped resampling and just restyled current map",
      ),
      "Go ahead and resample, I'll make a cup of coffee...",
      lambda: resample_active_map_for_em_half_angstrom(force=force, allow_large_output=True),
    )
    return None
  new_map_id = sharpen_blur_map_with_resampling(map_id, 0.0, resample_factor)
  if new_map_id == -1:
    return _restyle_active_map_without_resampling(
      map_id,
      contour_level_sigma,
      "Resampling failed",
    )
  style_resampled_em_map(new_map_id)
  set_contour_level_in_sigma(new_map_id, contour_level_sigma)
  set_map_displayed(new_map_id, 1)
  set_scroll_wheel_map(new_map_id)
  set_scrollable_map(new_map_id)
  if old_map_was_refinement_map:
    try:
      set_imol_refinement_map(new_map_id)
    except Exception:
      pass
  close_molecule(map_id)
  add_status_bar_text("Created EM-style 0.5 A/pixel map")
  return new_map_id


def make_masked_map_em(skip_non_em_confirmation=False):
  """Mask the active map by the active molecule and display the result in magenta."""
  map_id = _scrollable_map_or_status()
  if map_id is None:
    return None

  if not skip_non_em_confirmation and not _map_looks_em_like(map_id):
    generic_confirm_dialog(
      "Make masked map (EM)",
      "This is not an EM map. Are you sure?",
      "Cancel",
      lambda: None,
      "Continue",
      lambda: make_masked_map_em(skip_non_em_confirmation=True),
    )
    return None

  mol_id = _active_model_molecule_for_map_tools("masked map")
  if mol_id is None:
    return None

  current_contour_absolute = abs(get_contour_level_absolute(map_id))

  previous_mask_radius = None
  try:
    previous_mask_radius = map_mask_atom_radius()
  except Exception:
    previous_mask_radius = None

  try:
    set_map_mask_atom_radius(EM_MASKED_MAP_RADIUS)
    new_map_id = _mask_map_by_all_atoms(map_id, mol_id, 0)
  finally:
    if previous_mask_radius is not None:
      set_map_mask_atom_radius(previous_mask_radius)

  if new_map_id == -1 or new_map_id not in map_molecule_list():
    info_dialog("Failed to create masked map.")
    return None

  set_map_is_difference_map(new_map_id, 0)
  set_map_colour(
    new_map_id,
    EM_MASKED_MAP_COLOUR[0],
    EM_MASKED_MAP_COLOUR[1],
    EM_MASKED_MAP_COLOUR[2],
  )
  if current_contour_absolute > 0.0:
    set_contour_level_absolute(new_map_id, current_contour_absolute)
  set_map_displayed(new_map_id, 1)
  add_status_bar_text(
    "Created masked map #{0} from map #{1} around molecule #{2}".format(
      new_map_id,
      map_id,
      mol_id,
    )
  )
  return new_map_id


def make_pseudo_diff_map_em():
  """Backward-compatible wrapper for the older pseudo-diff menu action name."""
  return make_masked_map_em()


def _map_global_extent_radius(map_id):
  cell = map_cell(map_id)
  if not cell:
    return 20.0
  try:
    return max(cell[0], cell[1], cell[2]) / 2.0
  except Exception:
    return cell[0] / 2.0


def _save_global_view_settings(map_id):
  MAP_GLOBAL_VIEW_SETTINGS[map_id] = {
    "map_radius": _current_map_radius(map_id),
    "clipping_front": get_clipping_plane_front(),
    "clipping_back": get_clipping_plane_back(),
    "map_colour": list(map_colour_components_py(map_id)),
    "opacity": get_solid_density_surface_opacity(map_id),
  }


def _restore_global_view_settings(map_id, fallback_radius):
  saved = MAP_GLOBAL_VIEW_SETTINGS.get(map_id)
  if not saved:
    _set_map_radius_both(fallback_radius, map_id)
    return
  _set_map_radius_both(saved["map_radius"], map_id)
  set_clipping_front(saved["clipping_front"])
  set_clipping_back(saved["clipping_back"])


def _set_global_view_extent(map_id):
  current_radius = max(_current_map_radius(map_id), 1.0)
  global_radius = max(_map_global_extent_radius(map_id), current_radius)
  front = get_clipping_plane_front()
  back = get_clipping_plane_back()
  scale = global_radius / current_radius
  _set_map_radius_both(global_radius, map_id)
  set_clipping_front(front * scale)
  set_clipping_back(back * scale)


def toggle_global_map_view():
  map_id = _scrollable_map_or_status()
  if map_id is None:
    return None
  default_radius=20.0
  if map_is_difference_map(map_id)!=0:
    add_status_bar_text("Global surface view is only for non-difference maps")
    return None
  if not MAP_SURFACE_DISPLAY_STATE.get(map_id, False):
    local_appearance = MAP_LOCAL_APPEARANCE_STATE.get(map_id)
    _save_global_view_settings(map_id)
    set_draw_solid_density_surface(map_id,1)
    set_draw_map_standard_lines(map_id,0)
    set_flat_shading_for_solid_density_surface(0)
    set_map_fresnel_settings(map_id,1,0.0,0.5,2.0)
    set_solid_density_surface_opacity(
      map_id,
      MAP_SURFACE_OPACITY_STATE.get(map_id, get_solid_density_surface_opacity(map_id))
    )
    _set_global_view_extent(map_id)
    MAP_SURFACE_DISPLAY_STATE[map_id] = True
  else:
    saved_settings = MAP_GLOBAL_VIEW_SETTINGS.get(map_id, {})
    saved_colour = saved_settings.get("map_colour")
    set_draw_solid_density_surface(map_id,0)
    set_draw_map_standard_lines(map_id,1)
    saved_global_opacity = get_solid_density_surface_opacity(map_id)
    MAP_SURFACE_OPACITY_STATE[map_id] = saved_global_opacity
    set_solid_density_surface_opacity(map_id, saved_settings.get("opacity", 1.0))
    set_map_material_specular(map_id, 0.0, 64.0)
    set_map_fresnel_settings(map_id, 0, 0.0, 0.5, 2.0)
    if saved_colour and len(saved_colour) >= 3:
      set_map_colour(map_id,saved_colour[0],saved_colour[1],saved_colour[2])
    _restore_global_view_settings(map_id, default_radius)
    MAP_SURFACE_DISPLAY_STATE[map_id] = False
    
  
#Go to next residue in current polymer chain.
def next_res():
  reference_residue = _navigation_reference_residue()
  if not reference_residue:
    return None
  if reference_residue["source"] == "chain_start":
    return _go_to_navigation_residue(
      reference_residue["mol_id"],
      reference_residue["chain_id"],
      reference_residue["resno"],
      reference_residue["ins_code"],
    )
  mol_id = reference_residue["mol_id"]
  ch_id = reference_residue["chain_id"]
  serial_number = reference_residue.get("serial_number")
  if serial_number is None:
    serial_number = _residue_serial_number(
      mol_id,
      ch_id,
      reference_residue["resno"],
      reference_residue["ins_code"],
    )
  if serial_number < 0:
    add_status_bar_text("Could not resolve the current residue in its chain")
    return None
  n_residues = chain_n_residues(ch_id, mol_id)
  for next_serial in range(serial_number + 1, n_residues):
    next_resno = seqnum_from_serial_number(mol_id, ch_id, next_serial)
    next_ins_code = insertion_code_from_serial_number(mol_id, ch_id, next_serial)
    if _residue_is_polymer(mol_id, ch_id, next_resno, next_ins_code):
      return _go_to_navigation_residue(mol_id, ch_id, next_resno, next_ins_code, next_serial)
  add_status_bar_text("Already at the C-terminus of the current polymer chain")
  return None

#Go to previous residue in current polymer chain.
def prev_res():
  reference_residue = _navigation_reference_residue()
  if not reference_residue:
    return None
  if reference_residue["source"] == "chain_start":
    return _go_to_navigation_residue(
      reference_residue["mol_id"],
      reference_residue["chain_id"],
      reference_residue["resno"],
      reference_residue["ins_code"],
    )
  mol_id = reference_residue["mol_id"]
  ch_id = reference_residue["chain_id"]
  serial_number = reference_residue.get("serial_number")
  if serial_number is None:
    serial_number = _residue_serial_number(
      mol_id,
      ch_id,
      reference_residue["resno"],
      reference_residue["ins_code"],
    )
  if serial_number < 0:
    add_status_bar_text("Could not resolve the current residue in its chain")
    return None
  for previous_serial in range(serial_number - 1, -1, -1):
    previous_resno = seqnum_from_serial_number(mol_id, ch_id, previous_serial)
    previous_ins_code = insertion_code_from_serial_number(mol_id, ch_id, previous_serial)
    if _residue_is_polymer(mol_id, ch_id, previous_resno, previous_ins_code):
      return _go_to_navigation_residue(mol_id, ch_id, previous_resno, previous_ins_code, previous_serial)
  add_status_bar_text("Already at the N-terminus of the current polymer chain")
  return None
  
def sequence_context():
  residue = _active_residue_or_status()
  if residue is None:
    return None

  mol_id = residue[0]
  ch_id = residue[1]
  resnum = residue[2]
  ins_code = residue[3]
  resname = residue_name(mol_id, ch_id, resnum, ins_code)

  def get_aa_code(context_resnum):
    context_resname = residue_name(mol_id, ch_id, context_resnum, ins_code)
    if not context_resname:
      return "-"
    if len(context_resname) == 1:
      return context_resname
    aa_code = three_letter_code2single_letter(context_resname)
    if context_resname != "ALA" and aa_code == "A":
      return "X"
    return aa_code

  current_res = get_aa_code(resnum)
  minus_context = ""
  for i in range(1,10):
    minus_context = str(get_aa_code(resnum - i)) + minus_context
  plus_context = ""
  for i in range(1,10):
    plus_context = plus_context + str(get_aa_code(resnum + i))
  final_string="Residue:  " + resname + "  " + str(resnum) +"  Sequence context: ..." + minus_context + "[" + current_res + "]" + plus_context + "..."
  info_dialog(final_string)
  return 1
  
#Toggle display of active map
map_disp_flag={0:0}
map_disp_flag_cycle=0
def toggle_map_display():
  global map_disp_flag
  global map_disp_flag_cycle
  if map_disp_flag_cycle==0:
    for map_id in map_molecule_list():
      disp_value=map_is_displayed(map_id)
      map_disp_flag[map_id]=disp_value
      if disp_value==1:
        set_map_displayed(map_id,0) #If any maps are displayed, undisplay them.
    map_disp_flag_cycle=1
  elif map_disp_flag_cycle==1:
    disp_counter=0 
    for map_id in map_molecule_list():
      if map_id not in map_disp_flag: #if the map wasn't present in the previous cycle, assign a disp_value for it
        disp_value=map_is_displayed(map_id)
        map_disp_flag[map_id]=disp_value
      if map_disp_flag[map_id]==1:
        set_map_displayed(map_id,1) #Redisplay any maps that were displayed on the previous cycle.
      disp_counter=disp_counter+map_disp_flag[map_id] #test
    if disp_counter==0: #If no maps were displayed in the prior cycle, display all maps.
      for map_id in map_molecule_list(): 
        set_map_displayed(map_id,1) 
    map_disp_flag_cycle=0

#Colour active segment
def colour_active_segment():
  context = _active_segment_context_or_status()
  if not context:
    return None
  mol_id = context["mol_id"]
  selected = set(context["segment_specs"])
  colour_list = [(list(spec), 34) for spec in context["segment_specs"]]
  blank_list = [(spec, 0) for spec in _all_residue_specs_for_colouring(mol_id) if tuple(spec) not in selected]
  _apply_user_defined_residue_colours(mol_id, blank_list, colour_list)


def _all_residue_specs_for_colouring(mol_id):
  if callable(getattr(coot, "python_representation_kk", None)):
    try:
      return [list(spec) for spec in _model_coordinate_snapshot(mol_id)]
    except ValueError:
      pass
  residue_specs=[]
  for residue_entry in all_residues_with_serial_numbers(mol_id) or []:
    if not isinstance(residue_entry, list) or len(residue_entry) < 4:
      continue
    residue_spec = residue_entry[1:]
    chain_id = residue_spec_to_chain_id(residue_spec)
    resno = residue_spec_to_res_no(residue_spec)
    ins_code = residue_spec_to_ins_code(residue_spec)
    if chain_id is False or resno is False or ins_code is False:
      continue
    residue_specs.append([chain_id, resno, ins_code])
  return residue_specs


# ============================================================================
# Experimental custom-colouring helpers
# ============================================================================
#
# The old base-molecule user-defined colouring path leaks into ordinary display
# modes in this Coot 1.1 build. The only same-molecule route that looks
# plausible is to create additional representations using the dedicated
# user-defined bond-colour mode, then immediately clear the molecule's stored
# colour rules. Ramachandran colouring is the current test case for that path.
#

def _active_molecule_or_status():
  residue = _active_residue_or_status()
  if not residue:
    return None
  return residue[0]


def yellowify_carbons_in_active_molecule():
  """Use Coot's live bond-colour rotation path to yellowify active-molecule carbons."""
  mol_id = _active_molecule_or_status()
  if mol_id is None:
    return None
  previous_c_only_flag = get_colour_map_rotation_on_read_pdb_c_only_flag()
  try:
    set_colour_map_rotation_on_read_pdb_c_only_flag(1)
    # In Coot's bond-colour rotation code, 21 degrees is treated as the
    # default/yellow reference rotation for standard carbon colouring.
    set_bond_colour_rotation_for_molecule(mol_id, 21.0)
  finally:
    set_colour_map_rotation_on_read_pdb_c_only_flag(previous_c_only_flag)
  add_status_bar_text("Yellowified carbons in the active molecule")
  return 1


def _active_polymer_molecule_for_colouring(action_name):
  residue=_active_residue_or_status()
  if not residue:
    return None
  if not _residue_is_polymer(residue[0], residue[1], residue[2], residue[3]):
    info_dialog(action_name+" requires the active residue to be polymer.")
    return None
  return residue[0]


def _iter_chain_serial_numbers(mol_id, ch_id):
  try:
    n_residues = int(chain_n_residues(ch_id, mol_id))
  except Exception:
    return
  for serial_number in range(max(0, n_residues)):
    yield serial_number


def _clear_custom_colour_additional_representations(mol_id, clear_user_colours=True):
  handles = CUSTOM_COLOUR_ADDITIONAL_REPRESENTATIONS.pop(mol_id, [])
  for handle in handles:
    try:
      delete_additional_representation(mol_id, handle)
    except Exception:
      pass
  if clear_user_colours:
    try:
      clear_user_defined_atom_colours(mol_id)
    except Exception:
      pass


def _direct_user_defined_blank_rows_for_molecule(mol_id, blank_colour=0):
  del mol_id
  # Coot's CID "//" selects the whole molecule. Use one native rule rather
  # than building thousands of per-residue selections for large models.
  return [("//", blank_colour)]


def _graphics_to_legacy_user_defined_representation(mol_id, colour_rows, ca_only=False):
  if _supports_direct_user_defined_colouring():
    representation = (graphics_to_user_defined_atom_colours_representation if ca_only
                      else graphics_to_user_defined_atom_colours_all_atoms_representation)
    representation(mol_id)
    cycle_rep_flag[mol_id] = _representation_sequence().index(representation)
    return
  for item in colour_rows or []:
    if not isinstance(item, (list, tuple)) or len(item) < 2:
      continue
    residue_spec = item[0]
    if not isinstance(residue_spec, (list, tuple)) or len(residue_spec) < 3:
      continue
    chain_id = residue_spec[0]
    res_no = residue_spec[1]
    ins_code = residue_spec[2] or ""
    residue_name_here = residue_name(mol_id, chain_id, res_no, ins_code)
    if str(residue_name_here or "").upper() in WATER_RESIDUE_NAMES:
      graphics_to_user_defined_atom_colours_all_atoms_representation(mol_id)
      return
  graphics_to_user_defined_atom_colours_representation(mol_id)


def _set_user_defined_atom_colour_by_residue_atoms_py(mol_id, residue_specs_colour_index_tuple_list_py):
  """Apply residue colours atom-by-atom, avoiding the fragile selection-CID parser."""
  ensure_user_defined_colour_table()
  atom_colour_list = []
  for item in residue_specs_colour_index_tuple_list_py or []:
    if not isinstance(item, (list, tuple)) or len(item) < 2:
      continue
    residue_spec = item[0]
    colour_index = item[1]
    if not isinstance(residue_spec, (list, tuple)) or len(residue_spec) < 3:
      continue
    chain_id = residue_spec[0]
    res_no = residue_spec[1]
    ins_code = residue_spec[2] or ""
    atom_info = residue_info_py(mol_id, chain_id, res_no, ins_code)
    if not isinstance(atom_info, list):
      continue
    for atom in atom_info:
      try:
        atom_name = atom[0][0]
        alt_conf = atom[0][1]
      except Exception:
        continue
      atom_spec = [mol_id, chain_id, res_no, ins_code, atom_name, alt_conf]
      atom_colour_list.append(
        (atom_spec, _legacy_user_colour_index_to_coot_index(colour_index))
      )
  if not atom_colour_list:
    return 0
  try:
    coot.set_user_defined_atom_colour_py(mol_id, atom_colour_list)
  except Exception:
    traceback.print_exc()
    return 0
  return len(atom_colour_list)


def _make_custom_colour_additional_representation(mol_id, residue_spec):
  if not isinstance(residue_spec, (list, tuple)) or len(residue_spec) < 3:
    return None
  chain_id = residue_spec[0]
  res_no = residue_spec[1]
  ins_code = residue_spec[2] or ""
  draw_hydrogens_flag = 0
  if "draw_hydrogens_state" in globals():
    try:
      draw_hydrogens_flag = draw_hydrogens_state(mol_id)
    except Exception:
      draw_hydrogens_flag = 0
  handle = additional_representation_by_attributes(
    mol_id,
    chain_id,
    res_no,
    res_no,
    ins_code,
    CUSTOM_COLOUR_ADDITIONAL_REPRESENTATION_TYPE,
    CUSTOM_COLOUR_USER_DEFINED_BONDS_BOX_TYPE,
    CUSTOM_COLOUR_ADDITIONAL_BOND_WIDTH,
    draw_hydrogens_flag,
  )
  if handle == -1:
    return None
  return handle


def clear_custom_colour_representations_for_active_molecule():
  mol_id = _active_molecule_or_status()
  if mol_id is None:
    return None
  if _supports_direct_user_defined_colouring():
    restore_representation = _last_standard_representation_function(mol_id)
    # On fixed direct-colouring builds, clearing the colour rules while the
    # user-defined all-atoms representation is still active can trigger a crash
    # during redraw. Move back to the last remembered standard representation
    # first, then tear down the user-defined state.
    _restore_last_standard_representation(mol_id, restore_representation)
    _clear_custom_colour_additional_representations(mol_id, clear_user_colours=False)
    blank_rows = _direct_user_defined_blank_rows_for_molecule(mol_id)
    if blank_rows:
      _set_user_defined_atom_colour_by_residue_rows_py(mol_id, blank_rows)
  else:
    _clear_custom_colour_additional_representations(mol_id)
  if _supports_direct_user_defined_colouring():
    # Rebuild the same standard representation once the neutral overwrite has
    # happened. Avoiding clear_user_defined_atom_colours() on the fixed direct
    # path prevents Coot from regenerating with internal -1 colour indices.
    _restore_last_standard_representation(mol_id, restore_representation)
  add_status_bar_text("Cleared custom colour representations")


def _compact_direct_colour_rows(rows):
  """Resolve last-write-wins residue rules against a single neutral baseline."""
  resolved = {}
  for item in rows:
    if not isinstance(item, (list, tuple)) or len(item) < 2:
      continue
    spec, colour = item[:2]
    if not isinstance(spec, (list, tuple)) or len(spec) != 3:
      # Arbitrary CIDs may overlap residue rules; preserve their ordering.
      return rows
    resolved[tuple(spec)] = colour
  neutral = _legacy_user_colour_index_to_coot_index(0)
  return [(list(spec), colour) for spec, colour in resolved.items()
          if _legacy_user_colour_index_to_coot_index(colour) != neutral]


def _consume_colour_steps(iterator):
  while True:
    try:
      next(iterator)
    except StopIteration as result:
      return result.value


def _iter_direct_user_defined_residue_colours(mol_id, blank_res_list, colour_list, info_message=None):
  full_colour_list = [item[:2] for item in list(blank_res_list or []) + list(colour_list or [])
                     if isinstance(item, (list, tuple)) and len(item) >= 2]
  if not full_colour_list:
    if info_message:
      info_dialog(info_message)
    return None
  epoch = _MODEL_EDIT_EPOCHS.get(mol_id, 0)
  rows = _compact_direct_colour_rows(full_colour_list)
  yield "Preparing colours", 0.05, True
  snapshot = None
  # Explicit atom addressing is faster than thousands of full-model CID scans.
  # It addresses MODEL 1 only, so ensemble/ambiguous models retain the CID path.
  if (len(rows) >= 128 and all(isinstance(spec, (list, tuple)) for spec, _ in rows)
      and callable(getattr(coot, "python_representation_kk", None))
      and callable(getattr(coot, "set_user_defined_atom_colour_py", None))):
    try:
      snapshot = _model_coordinate_snapshot(mol_id)
    except ValueError:
      pass
  atom_rows = []
  if snapshot is not None:
    for index, (spec, colour) in enumerate(rows):
      record = snapshot.get(tuple(spec))
      if record is not None:
        native_colour = _legacy_user_colour_index_to_coot_index(colour)
        for atom in record["atoms"]:
          name, alt = atom[0][:2]
          atom_rows.append(([mol_id, *spec, name, alt], native_colour))
      if (index + 1) % 200 == 0:
        yield f"Preparing residues: {index + 1}/{len(rows)}", 0.05 + 0.45*(index+1)/len(rows), True
  yield "Applying colours", 0.5, True
  if not valid_model_molecule_qm(mol_id) or _MODEL_EDIT_EPOCHS.get(mol_id, 0) != epoch:
    raise RuntimeError("The model changed while preparing colours; please try again.")
  if snapshot is not None and snapshot != _model_coordinate_snapshot(mol_id):
    raise RuntimeError("The model changed while preparing colours; please try again.")
  # Once assignments start, complete the operation rather than leaving partial colours.
  yield "Applying colours", 0.5, False
  ensure_user_defined_colour_table()
  _clear_custom_colour_additional_representations(mol_id, clear_user_colours=False)
  neutral_rows = _direct_user_defined_blank_rows_for_molecule(mol_id)
  if not _set_user_defined_atom_colour_by_residue_rows_py(mol_id, neutral_rows, palette_ready=True):
    raise RuntimeError("Failed to initialise user-defined atom colours.")
  if snapshot is not None:
    for start in range(0, len(atom_rows), 4096):
      coot.set_user_defined_atom_colour_py(mol_id, atom_rows[start:start+4096])
      count = min(start+4096, len(atom_rows))
      yield f"Colouring atoms: {count}/{len(atom_rows)}", 0.5 + 0.45*count/len(atom_rows), False
  else:
    for start in range(0, len(rows), 64):
      if not _set_user_defined_atom_colour_by_residue_rows_py(mol_id, rows[start:start+64], palette_ready=True):
        raise RuntimeError("Failed to set user-defined atom colours.")
      count = min(start+64, len(rows))
      yield f"Colouring selections: {count}/{len(rows)}", 0.5 + 0.45*count/len(rows), False
  yield "Rebuilding display", 0.98, False
  _graphics_to_legacy_user_defined_representation(mol_id, neutral_rows + rows)
  yield "Colouring complete", 1.0, False
  if info_message:
    info_dialog(info_message)
  return len(rows) + 1


def _apply_direct_user_defined_residue_colours(mol_id, blank_res_list, colour_list, info_message=None):
  return _consume_colour_steps(_iter_direct_user_defined_residue_colours(
    mol_id, blank_res_list, colour_list, info_message))


def _iter_overlay_user_defined_residue_colours(mol_id, blank_res_list, colour_list, info_message=None):
  del blank_res_list
  if not colour_list:
    if info_message:
      info_dialog(info_message)
    return None
  yield "Assigning compatible atom colours", 0.05, False
  _clear_custom_colour_additional_representations(mol_id)
  clear_user_defined_atom_colours(mol_id)
  colour_count = _set_user_defined_atom_colour_by_residue_atoms_py(mol_id, colour_list)
  if not colour_count:
    _clear_custom_colour_additional_representations(mol_id)
    info_dialog("Failed to set user-defined atom colours for the requested residues.")
    return None
  handles = []
  CUSTOM_COLOUR_ADDITIONAL_REPRESENTATIONS[mol_id] = handles
  seen_residue_specs = set()
  yield "Building compatible representations", 0.2, False
  for index, (residue_spec, _colour_index) in enumerate(colour_list):
    residue_key = tuple(residue_spec[:3])
    if residue_key in seen_residue_specs:
      continue
    seen_residue_specs.add(residue_key)
    handle = _make_custom_colour_additional_representation(mol_id, residue_spec)
    if handle is not None:
      handles.append(handle)
    if (index + 1) % 4 == 0:
      yield f"Building representations: {index+1}/{len(colour_list)}", 0.2 + 0.75*(index+1)/len(colour_list), False
  clear_user_defined_atom_colours(mol_id)
  if not handles:
    info_dialog("Failed to create a custom-colour additional representation.")
    return None
  yield "Colouring complete", 1.0, False
  if info_message:
    info_dialog(info_message)
  return handles


def _apply_overlay_user_defined_residue_colours(mol_id, blank_res_list, colour_list, info_message=None):
  return _consume_colour_steps(_iter_overlay_user_defined_residue_colours(
    mol_id, blank_res_list, colour_list, info_message))


def _apply_user_defined_residue_colours(mol_id, blank_res_list, colour_list, info_message=None, defer_apply=False):
  if not defer_apply and _USER_DEFINED_COLOUR_JOB_ACTIVE:
    add_status_bar_text("Please wait for the current colouring operation to finish.")
    return None
  if _supports_direct_user_defined_colouring():
    if defer_apply:
      return _iter_direct_user_defined_residue_colours(mol_id, blank_res_list, colour_list, info_message)
    return _apply_direct_user_defined_residue_colours(
      mol_id,
      blank_res_list,
      colour_list,
      info_message=info_message,
    )
  steps = _iter_overlay_user_defined_residue_colours(mol_id, blank_res_list, colour_list, info_message)
  return steps if defer_apply else _consume_colour_steps(steps)


_USER_DEFINED_COLOUR_JOB_ACTIVE = False


def _start_user_defined_colour_action(mol_id, action, title):
  """Keep native graphics calls on the GTK thread, yielding between batches."""
  global _USER_DEFINED_COLOUR_JOB_ACTIVE
  if _USER_DEFINED_COLOUR_JOB_ACTIVE:
    add_status_bar_text("A colouring operation is already running.")
    return None
  glib = _coot_gui_repository_module("GLib")
  if glib is None or Gtk is None:
    return action(mol_id)
  window = Gtk.Window(title=title)
  window.set_default_size(360, 110)
  box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=8)
  for side in ("top", "bottom", "start", "end"):
    getattr(box, "set_margin_"+side)(12)
  label = Gtk.Label(label="Preparing colouring...")
  bar = Gtk.ProgressBar()
  cancel = Gtk.Button(label="Cancel")
  for widget in (label, bar, cancel):
    box.append(widget)
  window.set_child(box)
  state = {"cancelled": False, "cancellable": True, "iterator": None, "finished": False,
           "applying": False}
  epoch = _MODEL_EDIT_EPOCHS.get(mol_id, 0)
  _USER_DEFINED_COLOUR_JOB_ACTIVE = True
  def cancel_job(*_args):
    if state["cancellable"]:
      state["cancelled"] = True
    return True
  cancel.connect("clicked", cancel_job)
  window.connect("close-request", cancel_job)
  def finish():
    global _USER_DEFINED_COLOUR_JOB_ACTIVE
    if state["finished"]:
      return
    state["finished"] = True
    if state["iterator"] is not None:
      state["iterator"].close()
    _USER_DEFINED_COLOUR_JOB_ACTIVE = False
    window.destroy()
  def step():
    if state["finished"]:
      return False
    if state["cancelled"]:
      finish()
      add_status_bar_text("Colouring cancelled; existing colours retained.")
      return False
    try:
      if not valid_model_molecule_qm(mol_id) or _MODEL_EDIT_EPOCHS.get(mol_id, 0) != epoch:
        raise RuntimeError("The model changed during colouring; please try again.")
      if state["iterator"] is None:
        state["cancellable"] = False
        cancel.set_sensitive(False)
        result = action(mol_id, defer_apply=True)
        if result is None:
          finish()
          return False
        state["iterator"] = iter(result)
      message, fraction, cancellable = next(state["iterator"])
      state["cancellable"] = cancellable
      state["applying"] = state["applying"] or not cancellable
      label.set_text(message)
      cancel.set_sensitive(cancellable)
      if fraction is None:
        bar.pulse()
      else:
        bar.set_fraction(fraction)
      if fraction == 1.0:
        # Dismiss before any colour-key info dialog is presented.
        window.destroy()
      return True
    except StopIteration:
      finish()
      add_status_bar_text("Colouring complete.")
    except Exception as error:
      traceback.print_exc()
      if state["applying"] and valid_model_molecule_qm(mol_id):
        try:
          if not _supports_direct_user_defined_colouring():
            _clear_custom_colour_additional_representations(mol_id)
          _restore_last_standard_representation(mol_id)
        except Exception:
          traceback.print_exc()
      finish()
      info_dialog(f"Colouring stopped: {error}")
    return False
  label.set_text("Calculating colours / scores...")
  bar.pulse()
  window.present()
  # A timed idle-sized batch allows the progress window to paint before native work.
  glib.timeout_add(30, step)
  return window


def color_by_rama_native(mol_id, defer_apply=False):
  rama_results=all_molecule_ramachandran_score(mol_id)
  if not isinstance(rama_results, list) or len(rama_results) < 6:
    info_dialog("Unable to obtain Ramachandran scores.")
    return None
  scored_residues=rama_results[5]
  blank_colour=0
  rama_allowed_colour=27
  rama_outlier_colour=31
  blank_res_list=[]
  rama_colour_list=[]
  for residue_spec in _all_residue_specs_for_colouring(mol_id):
    blank_res_list.append((residue_spec, blank_colour))
  for item in scored_residues:
    if not isinstance(item, list) or len(item) < 3:
      continue
    residue_spec = _coot_residue_spec_from_spec(item[1])
    rama_score=item[2]
    if residue_spec is None:
      continue
    if rama_score < 0.002:
      rama_colour_list.append((list(residue_spec), rama_outlier_colour))
    elif rama_score < 0.02:
      rama_colour_list.append((list(residue_spec), rama_allowed_colour))
  return _apply_user_defined_residue_colours(
    mol_id,
    blank_res_list,
    rama_colour_list,
    "Ramachandran coloring:\n\nRed = outlier (<0.2%)\n\nOrange = allowed/disfavored (<2%)",
    defer_apply=defer_apply,
  )


def color_by_rama_native_for_active_residue():
  mol_id=_active_polymer_molecule_for_colouring("Ramachandran coloring")
  if mol_id is None:
    return None
  return _start_user_defined_colour_action(mol_id, color_by_rama_native, "Ramachandran colouring")


def color_by_density_fit_native(mol_id, defer_apply=False):
  map_id=imol_refinement_map()
  if map_id==-1:
    info_dialog("You need a refinement map for density-fit coloring.")
    return None
  residue_specs=all_residues_sans_water(mol_id)
  if not residue_specs:
    info_dialog("No residues found for density-fit coloring.")
    return None
  correlation_results=map_to_model_correlation_per_residue(mol_id, residue_specs, 0, map_id)
  blank_colour=0
  blank_res_list=[]
  density_colour_list=[]
  for residue_spec in _all_residue_specs_for_colouring(mol_id):
    blank_res_list.append((residue_spec, blank_colour))
  for item in correlation_results:
    if not isinstance(item, (list, tuple)) or len(item) < 2:
      continue
    residue_spec = _coot_residue_spec_from_spec(item[0])
    score=item[1]
    if residue_spec is None:
      continue
    try:
      score=float(score)
    except Exception:
      continue
    if math.isnan(score):
      continue
    colour_index=_density_fit_score_to_colour_index(score)
    density_colour_list.append((list(residue_spec), colour_index))
  return _apply_user_defined_residue_colours(
    mol_id,
    blank_res_list,
    density_colour_list,
    "Active molecule colored by model/map correlation.\n\n"
    + _colour_band_legend_text("CC", 0.0, 1.0),
    defer_apply=defer_apply,
  )


def color_by_density_fit_native_for_active_residue():
  mol_id=_active_polymer_molecule_for_colouring("Density-fit coloring")
  if mol_id is None:
    return None
  return _start_user_defined_colour_action(mol_id, color_by_density_fit_native, "Density-fit colouring")


def _ncs_difference_scores_from_ghosts(snapshot, ghosts):
  """Mean matched-atom displacement after Coot's peer-to-master NCS transform."""
  chains = {}
  for spec, record in snapshot.items():
    chains.setdefault(spec[0], {})[spec[1:]] = record
  atom_cache = {}
  def atom_positions(spec, record):
    if spec not in atom_cache:
      positions = {}
      for atom in record["atoms"]:
        name, alt = atom[0]
        point = tuple(float(value) for value in atom[2])
        if len(point) == 3 and all(math.isfinite(value) for value in point):
          positions.setdefault((str(name).strip(), str(alt or "").strip()), point)
      atom_cache[spec] = positions
    return atom_cache[spec]
  scores = {}
  seen_pairs = set()
  for ghost in ghosts:
    if not isinstance(ghost, (list, tuple)) or len(ghost) < 4:
      raise ValueError("Coot returned an unrecognised NCS operator.")
    peer, master, operator = ghost[1:4]
    if peer == master or (peer, master) in seen_pairs:
      continue
    if not isinstance(operator, (list, tuple)) or len(operator) != 2:
      raise ValueError("Coot could not calculate an NCS operator.")
    rotation, translation = (tuple(float(v) for v in part) for part in operator)
    if len(rotation) != 9 or len(translation) != 3 or not all(
        math.isfinite(value) for value in rotation + translation):
      raise ValueError("Coot returned an invalid NCS operator.")
    seen_pairs.add((peer, master))
    target_chain = chains.get(master, {})
    for key, record in chains.get(peer, {}).items():
      target = target_chain.get(key)
      name = str(record["name"]).strip().upper()
      if target is None or name in {"HOH", "WAT", "DOD"} or name != str(target["name"]).strip().upper():
        continue
      peer_spec, target_spec = (peer, *key), (master, *key)
      positions = atom_positions(peer_spec, record)
      target_positions = atom_positions(target_spec, target)
      total, count = 0.0, 0
      for atom_key, point in positions.items():
        other = target_positions.get(atom_key)
        if other is None:
          continue
        transformed = [sum(rotation[3*axis+j]*point[j] for j in range(3)) + translation[axis]
                       for axis in range(3)]
        total += math.sqrt(sum((a-b)**2 for a, b in zip(transformed, other)))
        count += 1
      if count:
        mean = total/count
        for spec in (peer_spec, target_spec):
          scores[spec] = max(scores.get(spec, 0.0), mean)
  return scores


def color_by_ncs_difference(mol_id, defer_apply=False):
  if mol_id not in model_molecule_list():
    info_dialog("You need an active model for NCS-difference coloring.")
    return None
  try:
    # ncs_chain_differences_py still slices the old residue-spec layout in
    # recent builds, losing the residue numbers. Read operators instead.
    ghosts = coot.ncs_ghosts_py(mol_id) or []
    if ghosts and any(not ghost[3] for ghost in ghosts):
      coot.make_ncs_ghosts_maybe(mol_id)
      ghosts = coot.ncs_ghosts_py(mol_id) or []
    if not ghosts:
      info_dialog("No NCS-related chain pairs were found for the active molecule.")
      return None
    snapshot = _model_coordinate_snapshot(mol_id)
    ncs_scores = _ncs_difference_scores_from_ghosts(snapshot, ghosts)
  except Exception as error:
    print("NCS-difference coloring failed:")
    traceback.print_exc()
    info_dialog(f"Could not calculate NCS differences: {error}")
    return None
  if not ncs_scores:
    info_dialog("NCS chains were found, but no matching residues/atoms could be compared.")
    return None

  blank_res_list = [(list(spec), 0) for spec in snapshot]

  ncs_colour_list=[]
  for residue_spec, mean_diff in ncs_scores.items():
    normalized_score=mean_diff/2.0
    if normalized_score < 0.0:
      normalized_score=0.0
    if normalized_score > 1.0:
      normalized_score=1.0
    colour_index=int(normalized_score*31+2)
    ncs_colour_list.append(([residue_spec[0], residue_spec[1], residue_spec[2]], colour_index))

  return _apply_user_defined_residue_colours(
    mol_id,
    blank_res_list,
    ncs_colour_list,
    "Active molecule colored by NCS difference.\n"
    "Mean matched-atom displacement after NCS superposition; each master residue shows its largest peer difference.\n"
    "Spectral scale: blue=0.0 A, red=2.0 A or greater; gray=not compared.",
    defer_apply=defer_apply,
  )


def color_by_ncs_difference_for_active_residue():
  mol_id=_active_polymer_molecule_for_colouring("NCS-difference coloring")
  if mol_id is None:
    return None
  return _start_user_defined_colour_action(mol_id, color_by_ncs_difference, "NCS-difference colouring")


def color_by_clash_score(mol_id, defer_apply=False):
  if mol_id not in model_molecule_list():
    info_dialog("You need an active model for clash coloring.")
    return None
  try:
    overlap_data=molecule_atom_overlaps(mol_id, -1)
  except TypeError:
    try:
      overlap_data=molecule_atom_overlaps(mol_id)
    except Exception:
      overlap_data=False
  except Exception:
    overlap_data=False
  if not isinstance(overlap_data, list):
    info_dialog("No clash data were available for the active molecule.")
    return None

  blank_colour=0
  blank_res_list=[]
  max_overlap_by_residue={}

  for residue_spec in _all_residue_specs_for_colouring(mol_id):
    blank_res_list.append((residue_spec, blank_colour))

  def accumulate_overlap_from_atom_spec(atom_spec, overlap_value):
    if not isinstance(atom_spec, list):
      return
    try:
      residue_spec=atom_spec_to_residue_spec(atom_spec)
    except Exception:
      return
    if not isinstance(residue_spec, list) or len(residue_spec) < 3:
      return
    residue_key=(residue_spec_to_chain_id(residue_spec),
                 residue_spec_to_res_no(residue_spec),
                 residue_spec_to_ins_code(residue_spec))
    previous_max=max_overlap_by_residue.get(residue_key, 0.0)
    if overlap_value > previous_max:
      max_overlap_by_residue[residue_key]=overlap_value

  for overlap_item in overlap_data:
    if not isinstance(overlap_item, dict):
      continue
    try:
      overlap_value=float(overlap_item.get('overlap-volume', 0.0))
    except Exception:
      continue
    atom_spec_1=overlap_item.get('atom-1-spec')
    atom_spec_2=overlap_item.get('atom-2-spec')
    accumulate_overlap_from_atom_spec(atom_spec_1, overlap_value)
    accumulate_overlap_from_atom_spec(atom_spec_2, overlap_value)

  if not max_overlap_by_residue:
    info_dialog("No clash data were available for the active molecule.")
    return None

  clash_colour_list=[]
  for residue_spec, max_overlap in max_overlap_by_residue.items():
    normalized_score=max_overlap/2.0
    if normalized_score < 0.0:
      normalized_score=0.0
    if normalized_score > 1.0:
      normalized_score=1.0
    colour_index=_density_fit_score_to_colour_index(1.0-normalized_score)
    clash_colour_list.append(([residue_spec[0], residue_spec[1], residue_spec[2]], colour_index))

  return _apply_user_defined_residue_colours(
    mol_id,
    blank_res_list,
    clash_colour_list,
    "Active molecule colored by per-residue maximum clash overlap.\n\n"
    + _colour_band_legend_text("overlap", 0.0, 2.0, reverse=True),
    defer_apply=defer_apply,
  )


def color_by_clash_score_for_active_molecule():
  mol_id=_active_molecule_or_status()
  if mol_id is None:
    return None
  return _start_user_defined_colour_action(mol_id, color_by_clash_score, "Clash colouring")

#Set refinement map to currently scrollable map
def set_map_to_scrollable_map():
  if map_molecule_list()!=[]:
    set_imol_refinement_map(scroll_wheel_map())
  else:
    info_dialog("You need a map for this to work.")
  
#Toggle display of active model
mol_disp_flag={0:0}
mol_disp_flag_cycle=0
def toggle_mol_display():
  global mol_disp_flag
  global mol_disp_flag_cycle
  if mol_disp_flag_cycle==0:
    for mol_id in model_molecule_list():
      mol_disp_value=mol_is_displayed(mol_id)
      mol_disp_flag[mol_id]=mol_disp_value
      if mol_disp_value==1:
        set_mol_displayed(mol_id,0)
    mol_disp_flag_cycle=1
  elif mol_disp_flag_cycle==1:
    disp_counter=0
    for mol_id in model_molecule_list():
      if mol_id not in mol_disp_flag:
        disp_value=mol_is_displayed(mol_id)
        mol_disp_flag[mol_id]=disp_value
      if mol_disp_flag[mol_id]==1:
        set_mol_displayed(mol_id,1)
      disp_counter=disp_counter+mol_disp_flag[mol_id]
    if disp_counter==0:
      for mol_id in model_molecule_list():
        set_mol_displayed(mol_id,1)
    mol_disp_flag_cycle=0

#Cycle representation mode forward/back
BASE_REPRESENTATION_SEQUENCE = [
  graphics_to_ca_plus_ligands_representation,
  graphics_to_ca_plus_ligands_and_sidechains_representation,
  graphics_to_ca_plus_ligands_sec_struct_representation,
  graphics_to_rainbow_representation,
  graphics_to_bonds_representation,
  graphics_to_b_factor_representation,
]
DEFAULT_REPRESENTATION_INDEX = 4
cycle_rep_flag = {}


def _representation_sequence():
  sequence = list(BASE_REPRESENTATION_SEQUENCE)
  if _supports_direct_user_defined_colouring():
    sequence.extend([
      graphics_to_user_defined_atom_colours_representation,
      graphics_to_user_defined_atom_colours_all_atoms_representation,
    ])
  return sequence

def _cycle_rep_flag_or_default(mol_id):
  representation_sequence = _representation_sequence()
  flag = cycle_rep_flag.get(mol_id, DEFAULT_REPRESENTATION_INDEX)
  # Seed an unknown model from Coot, then advance our requested position.
  # Re-reading on every keypress can pin the cycle if native mode reporting
  # lags or a native representation transition is a no-op.
  get_bond_type = getattr(coot, "get_graphics_molecule_bond_type", None)
  if mol_id not in cycle_rep_flag and callable(get_bond_type):
    bond_type_functions = {
      1: graphics_to_bonds_representation,
      15: graphics_to_bonds_representation,
      4: graphics_to_ca_plus_ligands_representation,
      17: graphics_to_ca_plus_ligands_and_sidechains_representation,
      7: graphics_to_ca_plus_ligands_sec_struct_representation,
      9: graphics_to_rainbow_representation,
      10: graphics_to_b_factor_representation,
      12: graphics_to_user_defined_atom_colours_all_atoms_representation,
      13: graphics_to_user_defined_atom_colours_representation,
    }
    try:
      current_function = bond_type_functions.get(get_bond_type(mol_id))
      if current_function in representation_sequence:
        flag = representation_sequence.index(current_function)
    except Exception:
      pass
  if flag not in range(len(representation_sequence)):
    flag = DEFAULT_REPRESENTATION_INDEX
  cycle_rep_flag[mol_id] = flag
  return flag

def _apply_representation_index(mol_id, flag):
  _representation_sequence()[flag](mol_id)
  cycle_rep_flag[mol_id] = flag
  return flag

def cycle_rep_up(mol_id,flag=None):
  representation_sequence = _representation_sequence()
  current_flag = _cycle_rep_flag_or_default(mol_id) if flag is None else flag
  next_flag = (current_flag + 1) % len(representation_sequence)
  return _apply_representation_index(mol_id, next_flag)


    
def cycle_rep_down(mol_id,flag=None):
  representation_sequence = _representation_sequence()
  current_flag = _cycle_rep_flag_or_default(mol_id) if flag is None else flag
  next_flag = (current_flag - 1) % len(representation_sequence)
  return _apply_representation_index(mol_id, next_flag)


def cycle_rep_up_current():
  residue = _active_residue_or_status()
  if not residue:
    return None
  mol_id = residue[0]
  cycle_rep_up(mol_id)


def cycle_rep_down_current():
  residue = _active_residue_or_status()
  if not residue:
    return None
  mol_id = residue[0]
  cycle_rep_down(mol_id)


# Embedded multimap tool: kept identical to multimap_refinement.py by tests.
# Load lazily in its own namespace so helper names cannot shadow trimmings.
_MULTIMAP_REFINEMENT_SOURCE = r'''"""Experimental, sequential multi-map refinement for Coot 1 / Python 3.

Load this file in Coot, then call multimap_refinement_gui(). Nothing runs on load.
Native edits are main-thread only; Pause/Stop take effect between refinement calls.
"""

import math
import traceback
from dataclasses import dataclass


_POLYMER_RESIDUE_NAMES = set("""A C G U T DA DC DG DT DU
PSU H2U 5MU OMU 4SU 1MA 2MA 6MA MIA 5MC OMC 1MG 2MG 7MG M2G OMG YG
ALA UNK ARG ASN ASP CYS GLU GLN GLY HIS ILE LEU LYS MET MSE PHE PRO SER THR TRP TYR VAL
P1L SEP TPO PTR CSO CME MLY MLZ HYP KCX ALY FME PYL SEC""".split())


def chain_residue_counts(records, polymer_specs=()):
    """Count polymer/non-polymer residues without confusing missing atoms with ligands."""
    polymer_specs = set(polymer_specs)
    counts = {}
    for spec, record in records.items():
        row = counts.setdefault(spec[0], [0, 0])
        polymer = spec in polymer_specs or record["name"].strip().upper() in _POLYMER_RESIDUE_NAMES
        row[0 if polymer else 1] += 1
    return counts


@dataclass(frozen=True)
class RefinementStage:
    map_id: int
    weight: float
    repeats: int = 1
    autofit_rotamers: bool = False

    def __post_init__(self):
        if type(self.map_id) is not int or self.map_id < 0:
            raise ValueError("Map must be a non-negative molecule number.")
        if not math.isfinite(self.weight) or self.weight <= 0:
            raise ValueError("Map weight must be finite and positive.")
        if type(self.repeats) is not int or self.repeats < 1:
            raise ValueError("Stage repeats must be a positive integer.")
        if type(self.autofit_rotamers) is not bool:
            raise ValueError("Autofit rotamers must be True or False.")


@dataclass(frozen=True)
class RefinementOptions:
    window_size: int = 11
    step_size: int = 5
    contact_radius: float = 4.0
    neighbour_padding: int = 1
    max_moving_residues: int = 0
    hydrogen_bonds: bool = False
    hbond_target: float = 2.9
    hbond_sigma: float = 0.15
    auto_pause: bool = False
    minimum_mean_cc: float = 0.2
    maximum_cc_drop: float = 0.1
    maximum_geometry: float = 5.0
    order: str = "window"
    follow_refinement: bool = True

    def __post_init__(self):
        for name in ("window_size", "step_size"):
            if type(getattr(self, name)) is not int or getattr(self, name) < 1:
                raise ValueError(name + " must be a positive integer.")
        if type(self.max_moving_residues) is not int or self.max_moving_residues < 0:
            raise ValueError("Moving residue limit must be a non-negative integer (0 = automatic).")
        if self.step_size > self.window_size:
            raise ValueError("Step size cannot exceed window size (would leave gaps).")
        if type(self.neighbour_padding) is not int or self.neighbour_padding < 0:
            raise ValueError("Neighbour padding must be a non-negative integer.")
        for name in ("contact_radius", "maximum_cc_drop"):
            if not math.isfinite(getattr(self, name)) or getattr(self, name) < 0:
                raise ValueError(name + " must be finite and non-negative.")
        for name in ("hbond_target", "hbond_sigma", "maximum_geometry"):
            if not math.isfinite(getattr(self, name)) or getattr(self, name) <= 0:
                raise ValueError(name + " must be finite and positive.")
        if not math.isfinite(self.minimum_mean_cc) or not -1 <= self.minimum_mean_cc <= 1:
            raise ValueError("Minimum CC must be between -1 and 1.")
        if self.order not in ("window", "stage"):
            raise ValueError("Order must be 'window' or 'stage'.")

    @property
    def moving_residue_limit(self):
        return self.max_moving_residues or max(500, 5 * self.window_size)


def _residue_spec(value):
    if not isinstance(value, (list, tuple)):
        raise ValueError("Invalid Coot residue specification.")
    if len(value) == 4:
        value = value[1:]
    if len(value) != 3 or not isinstance(value[0], str):
        raise ValueError("Invalid Coot residue specification: " + repr(value))
    return (value[0], int(value[1]), str(value[2] or ""))


def residue_windows(specs, window_size, step_size):
    """Use actual residue order, including insertion codes; do not leave a tail gap."""
    specs = tuple(_residue_spec(s) for s in specs)
    if not specs or len(set(specs)) != len(specs):
        raise ValueError("Range must contain distinct residues.")
    if len({s[0] for s in specs}) != 1 or list(specs) != sorted(specs):
        raise ValueError("Range must be ordered and within one chain.")
    if window_size < 1 or not 1 <= step_size <= window_size:
        raise ValueError("Invalid window/step sizes.")
    if len(specs) <= window_size:
        return [specs]
    starts = list(range(0, len(specs) - window_size + 1, step_size))
    if starts[-1] != len(specs) - window_size:
        starts.append(len(specs) - window_size)
    return [specs[i:i + window_size] for i in starts]


def refinement_windows(specs, window_size, step_size):
    """Split the scope by chain before constructing sliding windows."""
    specs = tuple(_residue_spec(s) for s in specs)
    if not specs or len(set(specs)) != len(specs):
        raise ValueError("Select at least one distinct polymer residue.")
    chains = {}
    for spec in specs:
        chains.setdefault(spec[0], []).append(spec)
    return [window for chain in chains.values()
            for window in residue_windows(chain, window_size, step_size)]


def _atom_key(spec):
    # Native atom specs may have a leading success/model field.
    if len(spec) == 6:
        spec = spec[1:]
    if len(spec) != 5:
        raise ValueError("Invalid atom specification: " + repr(spec))
    return (str(spec[0]), int(spec[1]), str(spec[2] or ""),
            str(spec[3]).strip(), str(spec[4] or ""))


def _distance(a, b):
    return math.sqrt(sum((x - y)**2 for x, y in zip(a, b)))


def _angle(a, b, c):
    u, v = [x-y for x, y in zip(a, b)], [x-y for x, y in zip(c, b)]
    length = math.sqrt(sum(x*x for x in u) * sum(x*x for x in v))
    if not length:
        return 0.0
    return math.degrees(math.acos(max(-1, min(1, sum(x*y for x, y in zip(u, v))/length))))


_CHI_ATOMS = {
    "ARG": ("N CA CB CG", "CA CB CG CD", "CB CG CD NE", "CG CD NE CZ"),
    "ASN": ("N CA CB CG", "CA CB CG OD1"),
    "ASP": ("N CA CB CG", "CA CB CG OD1"),
    "CYS": ("N CA CB SG",),
    "GLN": ("N CA CB CG", "CA CB CG CD", "CB CG CD OE1"),
    "GLU": ("N CA CB CG", "CA CB CG CD", "CB CG CD OE1"),
    "HIS": ("N CA CB CG", "CA CB CG ND1"),
    "ILE": ("N CA CB CG1", "CA CB CG1 CD1"),
    "LEU": ("N CA CB CG", "CA CB CG CD1"),
    "LYS": ("N CA CB CG", "CA CB CG CD", "CB CG CD CE", "CG CD CE NZ"),
    "MET": ("N CA CB CG", "CA CB CG SD", "CB CG SD CE"),
    "MSE": ("N CA CB CG", "CA CB CG SE", "CB CG SE CE"),
    "PHE": ("N CA CB CG", "CA CB CG CD1"),
    "PRO": ("N CA CB CG", "CA CB CG CD"),
    "SER": ("N CA CB OG",),
    "THR": ("N CA CB OG1",),
    "TRP": ("N CA CB CG", "CA CB CG CD1"),
    "TYR": ("N CA CB CG", "CA CB CG CD1"),
    "VAL": ("N CA CB CG1",),
}
_CHANGE_METRICS = {
    "ca": ("C-alpha RMSD", "A", (0.1, 0.25, 0.5, 1.0, 2.0)),
    "chi": ("Largest chi change", "deg", (5, 15, 30, 60, 120)),
    "backbone": ("Largest phi/psi/omega change", "deg", (5, 15, 30, 60, 120)),
}
_CHANGE_COLOURS = ((0.18, 0.36, 0.98), (0.08, 0.70, 0.90), (0.10, 0.78, 0.42),
                   (0.96, 0.90, 0.20), (0.98, 0.50, 0.06), (0.78, 0.00, 0.12))


def _comparison_coordinates(records):
    """Copy occupied coordinates only; no live native objects or molecule copies."""
    result = {}
    for spec, record in records.items():
        atoms = {}
        for atom in record["atoms"]:
            if float(atom[1][0]) <= 0:
                continue
            key = (str(atom[0][0]).strip(), str(atom[0][1] or ""))
            xyz = tuple(float(x) for x in atom[2])
            if key in atoms or len(xyz) != 3 or not all(math.isfinite(x) for x in xyz):
                raise ValueError("Ambiguous or invalid comparison coordinates: " + repr(spec))
            atoms[key] = xyz
        result[spec] = {"name": str(record["name"]).upper(), "atoms": atoms}
    return result


def _dihedral(points):
    def cross(a, b):
        return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
    bonds = [tuple(y-x for x, y in zip(a, b)) for a, b in zip(points, points[1:])]
    n1, n2 = cross(bonds[0], bonds[1]), cross(bonds[1], bonds[2])
    length = math.sqrt(sum(x*x for x in bonds[1]))
    if length < 1e-10 or min(sum(x*x for x in n1), sum(x*x for x in n2)) < 1e-20:
        return None
    axis = tuple(x/length for x in bonds[1])
    return math.degrees(math.atan2(sum(x*y for x, y in zip(cross(n1, n2), axis)),
                                   sum(x*y for x, y in zip(n1, n2))))


def _matched_points(initial, current, addresses):
    """Compare the same alternate identifier, allowing shared blank-alt atoms."""
    records = [model.get(spec) for model in (initial, current) for spec, _ in addresses]
    if any(record is None for record in records):
        return
    alts = {alt for record in records for _, alt in record["atoms"]} | {""}
    seen = set()
    for alt in sorted(alts):
        pairs = []
        for model in (initial, current):
            points = [model[spec]["atoms"].get((name, alt), model[spec]["atoms"].get((name, "")))
                      for spec, name in addresses]
            if any(point is None for point in points):
                break
            pairs.append(points)
        if len(pairs) == 2:
            key = tuple(tuple(points) for points in pairs)
            if key not in seen:
                seen.add(key)
                yield pairs


def _connected_pair(initial, current, left, right, backbone=False):
    if left[0] != right[0] or not (right[1] == left[1]+1 or
                                  (right[1] == left[1] and right[2] > left[2])):
        return False
    names, limits = (("C", "N"), (1.0, 2.0)) if backbone else (("CA", "CA"), (2.5, 4.5))
    return any(all(limits[0] <= _distance(*points) <= limits[1] for points in pair)
               for pair in _matched_points(initial, current, [(left, names[0]), (right, names[1])]))


def model_change_scores(initial, current, metric="ca", window_size=5):
    """Per-residue map-frame RMSD or maximum wrapped angular change.

    Inputs are detached records from _comparison_coordinates. Missing values
    are None. CA windows are centred, truncated at termini and chain breaks.
    Chi/backbone values take the largest change over comparable torsions/alts.
    """
    if metric not in _CHANGE_METRICS:
        raise ValueError("Unknown change metric: " + str(metric))
    if type(window_size) is not int or window_size < 1 or window_size % 2 != 1:
        raise ValueError("C-alpha window size must be a positive odd integer.")
    specs = sorted(s for s in current if {"N", "CA"} <= {name for name, _ in current[s]["atoms"]})
    values = {s: None for s in specs}
    if metric == "ca":
        segments = []
        for spec in specs:
            if not segments or not _connected_pair(initial, current, segments[-1][-1], spec):
                segments.append([])
            segments[-1].append(spec)
        for segment in segments:
            squares = []
            for spec in segment:
                pairs = list(_matched_points(initial, current, [(spec, "CA")]))
                squares.append(sum(_distance(a[0], b[0])**2 for a, b in pairs)/len(pairs) if pairs else None)
            sums, counts = [0.0], [0]
            for value in squares:
                sums.append(sums[-1] + (value or 0.0))
                counts.append(counts[-1] + (value is not None))
            half = window_size//2
            for i, spec in enumerate(segment):
                lo, hi = max(0, i-half), min(len(segment), i+half+1)
                count = counts[hi]-counts[lo]
                if squares[i] is not None and count:
                    values[spec] = math.sqrt(max(0, sums[hi]-sums[lo])/count)
        return values
    for i, spec in enumerate(specs):
        torsions = []
        if metric == "chi":
            if spec not in initial or initial[spec]["name"] != current[spec]["name"]:
                continue
            name = current[spec]["name"]
            for j, atoms in enumerate(_CHI_ATOMS.get(name, ())):
                period = 180 if (name, j) in (("ASP", 1), ("GLU", 2), ("PHE", 1), ("TYR", 1)) else 360
                torsions.append(([(spec, atom) for atom in atoms.split()], period))
        else:
            if i and _connected_pair(initial, current, specs[i-1], spec, backbone=True):
                torsions.append(([(specs[i-1], "C"), (spec, "N"), (spec, "CA"), (spec, "C")], 360))
            if i+1 < len(specs) and _connected_pair(initial, current, spec, specs[i+1], backbone=True):
                other = specs[i+1]
                torsions.extend((([(spec, "N"), (spec, "CA"), (spec, "C"), (other, "N")], 360),
                                 ([(spec, "CA"), (spec, "C"), (other, "N"), (other, "CA")], 360)))
        deltas = []
        for addresses, period in torsions:
            for pair in _matched_points(initial, current, addresses):
                if metric == "backbone" and any(
                        addresses[k][0] != addresses[k+1][0] and not 1.0 <= _distance(points[k], points[k+1]) <= 2.0
                        for points in pair for k in range(3)):
                    continue
                a, b = (_dihedral(points) for points in pair)
                if a is not None and b is not None:
                    deltas.append(abs((b-a+period/2) % period-period/2))
        if deltas:
            values[spec] = max(deltas)
    return values


def change_colour_legend(metric):
    title, unit, limits = _CHANGE_METRICS[metric]
    bands = (0,) + limits
    lines = [title + " (" + unit + ")"]
    for i, name in enumerate(("Blue", "Cyan", "Green", "Yellow", "Orange", "Red")):
        lines.append("%s: %g-%g" % (name, bands[i], bands[i+1]) if i < 5 else "Red: >= %g" % bands[i])
    return "\n".join(lines + ["Grey: unprocessed or no comparable atoms/torsions"])


def _change_colour_table():
    # Align positional CA lookup with explicit atom-colour keys, without
    # overriding standard 0..59 colours or clearing native colours to -1.
    neutral = 120
    palette = [(neutral, [0.72, 0.72, 0.72]) for _ in range(neutral+1)]
    palette.extend((neutral+1+i, list(rgb)) for i, rgb in enumerate(_CHANGE_COLOURS))
    return palette


def backbone_hbond_pairs(records, moving):
    """Conservative O...N candidates, without changing global H-bond settings.

    Backbone only, no proline donors, 2.6--3.3 A, C-O-N and CA-N-O >=110 deg.
    Separate alternate conformers are never joined. At least one endpoint moves.
    """
    donors, acceptors = [], []
    for spec, record in records.items():
        atoms = record["atoms"]
        for (name, alt), (native, xyz) in atoms.items():
            if name == "N" and record["name"] != "PRO":
                ca = atoms.get(("CA", alt)) or atoms.get(("CA", ""))
                if ca:
                    donors.append((spec, alt, native, xyz, ca[1]))
            if name == "O":
                carbon = atoms.get(("C", alt)) or atoms.get(("C", ""))
                if carbon:
                    acceptors.append((spec, alt, native, xyz, carbon[1]))
    # Local spatial bins keep this cheap even for large contact shells.
    bins = {}
    for acceptor in acceptors:
        key = tuple(math.floor(x / 3.3) for x in acceptor[3])
        bins.setdefault(key, []).append(acceptor)
    result = []
    for ds, da, dn, dp, dca in donors:
        key = tuple(math.floor(x / 3.3) for x in dp)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for os, oa, on, op, oc in bins.get((key[0]+dx, key[1]+dy, key[2]+dz), []):
                        if ds == os or (ds not in moving and os not in moving):
                            continue
                        if da and oa and da != oa:
                            continue
                        if ds[0] == os[0] and abs(ds[1] - os[1]) <= 1:
                            continue
                        if (2.6 <= _distance(dp, op) <= 3.3 and
                                _angle(oc, op, dp) >= 110 and _angle(dca, dp, op) >= 110):
                            result.append((dn, on))
    return result


class CootRefinementBackend:
    """Coot adapter. Settings and only our own temporary restraints are scoped per step."""

    def __init__(self, api=None):
        if api is None:
            import coot
            api = coot
        self.api = api
        required = ("python_representation_kk", "refine_residues_py", "accept_regularizement",
                    "matrix_state", "set_matrix", "imol_refinement_map", "set_imol_refinement_map",
                    "refinement_immediate_replacement_state", "set_refinement_immediate_replacement",
                    "is_valid_model_molecule", "is_valid_map_molecule",
                    "residue_info_py", "residue_name")
        missing = [name for name in required if not callable(getattr(api, name, None))]
        if missing:
            raise RuntimeError("Unsupported Coot build; missing: " + ", ".join(missing))
        self.chain_specs = {}
        self.polymer_specs = set()
        self.positions = {}
        self._grid_radius = None
        self._grid = {}
        self._grid_keys = {}

    def snapshot(self, imol):
        if not self.api.is_valid_model_molecule(imol):
            raise ValueError("Model has been closed or is not a coordinate molecule.")
        count = getattr(self.api, "n_models", None)
        if callable(count) and count(imol) != 1:
            raise ValueError("A single-model coordinate molecule is required.")
        models = self.api.python_representation_kk(imol)
        if not models or len(models) != 1:
            raise ValueError("A single-model coordinate molecule is required.")
        records = {}
        self.chain_specs = {}
        self.polymer_specs = set()
        self.positions = {}
        self._grid_radius = None
        self._grid = {}
        self._grid_keys = {}
        for chain, residues in models[0]:
            for number, ins, name, atoms in residues:
                spec = (str(chain), int(number), str(ins or ""))
                if spec in records:
                    raise ValueError("Ambiguous residue identifiers in model.")
                records[spec] = {"name": name, "atoms": atoms}
                self.positions[spec] = self._positions(atoms)
                names = {atom[0][0].strip() for atom in atoms}
                if {"N", "CA", "C"} <= names or {"C3'", "C4'"} <= names:
                    self.polymer_specs.add(spec)
                self.chain_specs.setdefault(spec[0], []).append(spec)
        for specs in self.chain_specs.values():
            specs.sort()
        return records

    def colour_changes(self, imol, records, values, metric):
        graphics_state = getattr(self.api, "use_graphics_interface_state", None)
        if callable(graphics_state) and graphics_state() == 0:
            raise RuntimeError("Change colouring requires Coot's graphical interface.")
        ca_representation = "graphics_to_user_defined_atom_colours_representation"
        atom_representation = "graphics_to_user_defined_atom_colours_all_atoms_representation"
        representations = ((ca_representation, atom_representation) if metric == "chi"
                           else (atom_representation, ca_representation))
        for name in ("set_user_defined_colours_py", "set_user_defined_atom_colour_py", *representations):
            if not callable(getattr(self.api, name, None)):
                raise RuntimeError("This Coot build cannot display change colours; missing " + name)
        if not values:
            raise ValueError("No comparable atoms/torsions for this metric.")
        neutral = 120
        limits = _CHANGE_METRICS[metric][2]
        self.api.set_user_defined_colours_py(_change_colour_table())
        rows = []
        for spec, record in records.items():
            value = values.get(spec)
            colour = neutral if value is None else neutral+1+sum(value >= cutoff for cutoff in limits)
            for atom in record["atoms"]:
                rows.append(([imol, *spec, *atom[0][:2]], colour))
                if len(rows) == 4096:
                    self.api.set_user_defined_atom_colour_py(imol, rows)
                    rows = []
        if rows:
            self.api.set_user_defined_atom_colour_py(imol, rows)
        # Refresh both views from the new assignments, ending in the metric's
        # usual view so switching representations cannot retain an older metric.
        for representation in representations:
            getattr(self.api, representation)(imol)
        self.redraw()

    @staticmethod
    def _positions(atoms):
        points = []
        for atom in atoms:
            if float(atom[1][0]) > 0:
                xyz = tuple(float(x) for x in atom[2])
                if len(xyz) != 3 or not all(math.isfinite(x) for x in xyz):
                    raise ValueError("Invalid atom coordinates.")
                points.append(xyz)
        return points

    def _index_positions(self, spec):
        radius = self._grid_radius
        keys = set()
        for xyz in self.positions[spec]:
            key = tuple(math.floor(x / radius) for x in xyz)
            self._grid.setdefault(key, {}).setdefault(spec, []).append(xyz)
            keys.add(key)
        self._grid_keys[spec] = keys

    def _refresh_positions(self, imol, specs):
        for spec in specs:
            self.positions[spec] = self._positions(self.api.residue_info_py(imol, *spec) or [])
            if self._grid_radius is not None:
                for key in self._grid_keys.get(spec, ()):
                    bucket = self._grid[key]
                    bucket.pop(spec, None)
                    if not bucket:
                        del self._grid[key]
                self._index_positions(spec)

    def validate(self, imol, stages, options, specs):
        self._check_idle()
        records = self.snapshot(imol)
        if not self.api.is_valid_map_molecule(self.api.imol_refinement_map()):
            raise ValueError("Select a refinement map in Coot before starting (needed for state restoration).")
        if any(s not in records for s in specs):
            raise ValueError("A requested residue no longer exists.")
        for stage in stages:
            if not self.api.is_valid_map_molecule(stage.map_id):
                raise ValueError("Map #%s has been closed or is not a map." % stage.map_id)
        if options.auto_pause and not callable(getattr(self.api, "map_to_model_correlation_per_residue_py", None)):
            raise RuntimeError("CC guard is unavailable in this Coot build.")
        if options.hydrogen_bonds:
            for name in ("add_extra_bond_restraint", "delete_extra_restraint_py", "list_extra_restraints_py"):
                if not callable(getattr(self.api, name, None)):
                    raise RuntimeError("Temporary H-bond restraints require " + name)
        if any(stage.autofit_rotamers for stage in stages):
            for name in ("auto_fit_best_rotamer", "rotamer_auto_fit_do_post_refine_state",
                         "set_rotamer_auto_fit_do_post_refine"):
                if not callable(getattr(self.api, name, None)):
                    raise RuntimeError("Rotamer autofitting requires " + name)

    def _check_idle(self):
        ongoing = getattr(self.api, "refinement_already_ongoing_p", None)
        if callable(ongoing) and ongoing():
            raise RuntimeError("Finish/accept/reject the existing Coot refinement first.")

    def moving_selection(self, imol, seed, options):
        selected = set(seed)
        if options.contact_radius:
            # Native contact wrappers share borrowed Python references on some builds.
            # Use a spatial index and update only moved residues after each refinement.
            radius = options.contact_radius
            if self._grid_radius != radius:
                self._grid_radius, self._grid, self._grid_keys = radius, {}, {}
                for spec in self.positions:
                    self._index_positions(spec)
            radius_sq = radius * radius
            for spec in seed:
                for xyz in self.positions[spec]:
                    key = tuple(math.floor(x / radius) for x in xyz)
                    for dx in (-1, 0, 1):
                        for dy in (-1, 0, 1):
                            for dz in (-1, 0, 1):
                                for other, points in self._grid.get((key[0]+dx, key[1]+dy, key[2]+dz), {}).items():
                                    if other not in selected and any(sum((a-b)**2 for a, b in zip(xyz, p)) <= radius_sq for p in points):
                                        selected.add(other)
        # Pad contact residues by actual chain order, retaining waters/ligands as contacts.
        contacts = set(selected)
        for chain, ordered in self.chain_specs.items():
            ordered = [s for s in ordered if s in self.polymer_specs]
            indices = {s: i for i, s in enumerate(ordered)}
            for spec in contacts:
                if spec[0] == chain and spec in indices:
                    i = indices[spec]
                    selected.update(ordered[max(0, i-options.neighbour_padding):i+options.neighbour_padding+1])
        if len(selected) > options.moving_residue_limit:
            raise ValueError("Moving selection has %d residues (window %d + contacts %d; limit %d). "
                             "Increase the moving residue limit in Options, or reduce the contact radius/padding." %
                             (len(selected), len(seed), len(selected)-len(seed), options.moving_residue_limit))
        for spec in selected:
            if not self.api.residue_info_py(imol, *spec):
                raise ValueError("Moving residue disappeared: " + repr(spec))
        return sorted(selected)

    def local_records(self, imol, specs):
        records = {}
        for spec in specs:
            atoms = {}
            for atom in self.api.residue_info_py(imol, *spec) or []:
                name, alt = atom[0]
                xyz = tuple(float(x) for x in atom[2])
                if not all(math.isfinite(x) for x in xyz):
                    raise ValueError("Non-finite coordinates in " + repr(spec))
                if float(atom[1][0]) <= 0:
                    continue
                key = (name.strip(), str(alt or ""))
                if key in atoms:
                    raise ValueError("Ambiguous atom identity in " + repr(spec))
                atoms[key] = (list(spec) + [name, str(alt or "")], xyz)
            records[spec] = {"name": self.api.residue_name(imol, *spec), "atoms": atoms}
        return records

    def correlation(self, imol, specs, map_id):
        rows = self.api.map_to_model_correlation_per_residue_py(imol, [list(s) for s in specs], 0, map_id)
        scores = {}
        for spec, value in rows or []:
            spec = _residue_spec(spec)
            if spec in specs and math.isfinite(float(value)):
                scores[spec] = float(value)
        if len(scores) != len(specs):
            raise RuntimeError("CC guard could score only %d/%d residues." % (len(scores), len(specs)))
        return sum(scores.values()) / len(scores)

    def autofit_rotamers(self, imol, seed, stage):
        """Fit the window (not its contact shell), preserving search/backrub settings."""
        self._check_idle()
        api = self.api
        if not api.is_valid_model_molecule(imol) or not api.is_valid_map_molecule(stage.map_id):
            raise ValueError("Model or stage map has been closed.")
        previous_map = api.imol_refinement_map()
        post_refine = api.rotamer_auto_fit_do_post_refine_state()
        if not api.is_valid_map_molecule(previous_map):
            raise ValueError("Select a refinement map before starting.")
        report = {"attempted": 0, "failed": []}
        amino_acids = {"ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "HIS", "ILE", "LEU",
                      "LYS", "MET", "MSE", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}
        try:
            if api.set_imol_refinement_map(stage.map_id) != stage.map_id:
                raise RuntimeError("Could not set rotamer-fitting map.")
            api.set_rotamer_auto_fit_do_post_refine(0)
            for spec in seed:
                if api.residue_name(imol, *spec).strip().upper() not in amino_acids:
                    continue
                atoms = api.residue_info_py(imol, *spec) or []
                alts = sorted({str(a[0][1]) for a in atoms if a[0][1] and float(a[1][0]) > 0})
                if not any(float(a[1][0]) > 0 for a in atoms):
                    continue
                for alt in alts or [""]:
                    score = api.auto_fit_best_rotamer(imol, *spec, alt, stage.map_id, 1, 0.3)
                    report["attempted"] += 1
                    if not isinstance(score, (int, float)) or not math.isfinite(score) or score <= -99.89:
                        report["failed"].append(tuple(spec) + (alt,))
            return report
        finally:
            errors = []
            # Fitting/backrub may move neighbours or delete clashing waters.
            for function, args in ((api.set_rotamer_auto_fit_do_post_refine, (post_refine,)),
                                   (api.set_imol_refinement_map, (previous_map,)),
                                   (self.snapshot, (imol,)), (self.redraw, ())):
                try:
                    function(*args)
                except Exception as error:
                    errors.append(str(error))
            if errors:
                raise RuntimeError("Rotamer cleanup failed: " + "; ".join(errors))

    def navigate(self, imol, seed):
        """Centre on backbone anchors, without changing orientation or zoom."""
        anchors = []
        for spec in seed:
            atoms = self.api.residue_info_py(imol, *spec) or []
            backbone = [a for a in atoms if a[0][0].strip() in ("CA", "C4'") and float(a[1][0]) > 0]
            if backbone:
                anchors.append(tuple(float(x) for x in max(backbone, key=lambda a: float(a[1][0]))[2]))
            elif self.positions.get(spec):
                points = self.positions[spec]
                anchors.append(tuple(sum(p[i] for p in points) / len(points) for i in range(3)))
        if not anchors:
            raise ValueError("No occupied atoms remain in the working region.")
        centre = tuple(sum(p[i] for p in anchors) / len(anchors) for i in range(3))
        if not all(math.isfinite(x) for x in centre):
            raise ValueError("Invalid coordinates in the working region.")
        self.api.set_rotation_centre(*centre)
        self.redraw()

    def redraw(self):
        draw = getattr(self.api, "graphics_draw", None)
        if callable(draw):
            draw()

    def refine(self, imol, specs, stage, options):
        api = self.api
        self._check_idle()
        if not api.is_valid_model_molecule(imol) or not api.is_valid_map_molecule(stage.map_id):
            raise ValueError("Model or stage map has been closed.")
        previous = (api.imol_refinement_map(), api.matrix_state(), api.refinement_immediate_replacement_state())
        if not api.is_valid_map_molecule(previous[0]):
            raise ValueError("Select a refinement map before starting.")
        owned = []
        attempted = False
        try:
            if options.hydrogen_bonds:
                existing = {frozenset((_atom_key(row[1]), _atom_key(row[2])))
                            for row in api.list_extra_restraints_py(imol) or [] if row[0] == "bond"}
                for a, b in backbone_hbond_pairs(self.local_records(imol, specs), set(specs)):
                    key = frozenset((_atom_key(a), _atom_key(b)))
                    if key in existing:
                        continue
                    index = api.add_extra_bond_restraint(imol, *a, *b, options.hbond_target, options.hbond_sigma)
                    if index < 0:
                        raise RuntimeError("Could not add temporary H-bond restraint.")
                    owned.append(["bond", a, b])
                    existing.add(key)
            if api.set_imol_refinement_map(stage.map_id) != stage.map_id:
                raise RuntimeError("Could not set refinement map.")
            api.set_matrix(stage.weight)
            api.set_refinement_immediate_replacement(1)
            attempted = True
            result = api.refine_residues_py(imol, [list(s) for s in specs])
            if result is False or result is None:
                raise RuntimeError("Coot could not set up refinement (check dictionaries/selection).")
            ongoing = getattr(api, "refinement_already_ongoing_p", None)
            if not callable(ongoing) or ongoing():
                accept = getattr(api, "accept_moving_atoms_py", None)
                if callable(accept):
                    accepted = accept()
                    if accepted:
                        result = accepted
                else:
                    api.accept_regularizement()
            self._refresh_positions(imol, specs)
            self.redraw()
            return result, len(owned)
        except Exception:
            ongoing = getattr(api, "refinement_already_ongoing_p", None)
            clear = getattr(api, "clear_up_moving_atoms", None)
            if attempted and callable(ongoing) and ongoing() and callable(clear):
                clear()
            raise
        finally:
            # Run every cleanup even if another cleanup fails; surface errors to the runner.
            cleanup_errors = []
            actions = [(api.delete_extra_restraint_py, (imol, row)) for row in owned] if owned else []
            actions += [(api.set_refinement_immediate_replacement, (previous[2],)),
                        (api.set_matrix, (previous[1],)), (api.set_imol_refinement_map, (previous[0],))]
            for function, args in actions:
                try:
                    function(*args)
                except Exception as error:
                    cleanup_errors.append(str(error))
            if cleanup_errors:
                raise RuntimeError("Cleanup failed: " + "; ".join(cleanup_errors))


_MULTIMAP_RUNNING = None


class MultimapRefinementJob:
    """One bounded native refinement per tick; accepted edits remain on stop/error."""

    def __init__(self, imol, specs, stages, options=None, backend=None, on_update=None):
        self.imol = imol
        self.options = options or RefinementOptions()
        self.specs = tuple(_residue_spec(s) for s in specs)
        self.stages = tuple(stages)
        if not self.stages:
            raise ValueError("Add at least one refinement stage.")
        self.windows = refinement_windows(self.specs, self.options.window_size, self.options.step_size)
        self.backend = backend or CootRefinementBackend()
        self.backend.validate(imol, self.stages, self.options, self.specs)
        stage_steps = [(i, r) for i, s in enumerate(self.stages) for r in range(s.repeats)]
        if self.options.order == "window":
            self.tasks = [(w, i, r) for w in range(len(self.windows)) for i, r in stage_steps]
        else:
            self.tasks = [(w, i, r) for i, r in stage_steps for w in range(len(self.windows))]
        self.index = 0
        self.state = "ready"
        self.message = "Ready: %d windows, %d refinement steps." % (len(self.windows), len(self.tasks))
        self.history = []
        self.processed_residues = set()
        self.on_update = on_update
        self._busy = False
        self._requested = None
        self._prepared_index = None
        self.initial_coordinates = None

    def _notify(self):
        print("Multimap refinement: " + self.message)
        if self.on_update:
            self.on_update(self)

    def start(self):
        global _MULTIMAP_RUNNING
        if self.state not in ("ready", "paused"):
            return
        if self.index == len(self.tasks):
            self._settle("complete", "Complete; no remaining refinement steps.")
            return
        if _MULTIMAP_RUNNING not in (None, self):
            raise RuntimeError("Another multimap job is running; pause or stop it first.")
        self.backend.validate(self.imol, self.stages, self.options, self.specs)
        if self.initial_coordinates is None:
            snapshot = getattr(self.backend, "snapshot", None)
            if callable(snapshot):
                self.initial_coordinates = _comparison_coordinates(snapshot(self.imol))
        _MULTIMAP_RUNNING = self
        self.state, self._requested = "running", None
        self._prepared_index = None
        self.message = "Running %d/%d." % (self.index, len(self.tasks))
        self._notify()

    def visualize_changes(self, metric="ca", window_size=5):
        """Colour the current model against this run's initial coordinates."""
        if self.state not in ("paused", "complete", "stopped", "failed") or self.initial_coordinates is None:
            raise RuntimeError("Pause, finish or stop this run before displaying changes.")
        if _MULTIMAP_RUNNING is not None:
            raise RuntimeError("Stop the other multimap run before displaying changes.")
        self.backend._check_idle()
        records = self.backend.snapshot(self.imol)
        values = model_change_scores(self.initial_coordinates, _comparison_coordinates(records), metric, window_size)
        values = {spec: value if spec in self.processed_residues else None
                  for spec, value in values.items()}
        self.backend.colour_changes(self.imol, records, values, metric)
        return values

    def _settle(self, state, message):
        global _MULTIMAP_RUNNING
        self.state, self.message = state, message
        if _MULTIMAP_RUNNING is self:
            _MULTIMAP_RUNNING = None
        self._notify()

    def pause(self):
        if self.state == "running":
            if self._busy:
                self._requested = "paused"
            else:
                self._settle("paused", "Paused at %d/%d." % (self.index, len(self.tasks)))

    def stop(self):
        if self.state in ("ready", "running", "paused"):
            if self._busy:
                self._requested = "stopped"
            else:
                self._settle("stopped", "Stopped; completed edits retained. Use Coot Undo if needed.")

    def prepare_step(self):
        """Publish the working region before a separate callback runs native work."""
        if self.state != "running" or self._busy:
            return False
        if self._prepared_index == self.index:
            return True
        try:
            w, si, repeat = self.tasks[self.index]
            seed, stage = self.windows[w], self.stages[si]
            navigate = getattr(self.backend, "navigate", None)
            if self.options.follow_refinement and callable(navigate):
                navigate(self.imol, seed)
            self.message = "Refining %d/%d | %s:%d%s-%d%s | Stage %d/%d, map #%d, weight %g | Repeat %d/%d" % (
                self.index+1, len(self.tasks), seed[0][0], seed[0][1], seed[0][2], seed[-1][1], seed[-1][2],
                si+1, len(self.stages), stage.map_id, stage.weight, repeat+1, stage.repeats)
            if stage.autofit_rotamers:
                self.message += " | Autofit rotamers, then refine"
            self._prepared_index = self.index
            self._notify()
        except Exception as error:
            traceback.print_exc()
            self._settle("failed", "Failed: %s. Earlier edits remain; use Coot Undo if needed." % error)
        return self.state == "running"

    def tick(self):
        if not self.prepare_step():
            return False
        self._busy = True
        try:
            w, si, repeat = self.tasks[self.index]
            seed, stage = self.windows[w], self.stages[si]
            reference = self.stages[-1].map_id
            before = self.backend.correlation(self.imol, seed, reference) if self.options.auto_pause else None
            rotamers = self.backend.autofit_rotamers(self.imol, seed, stage) if stage.autofit_rotamers else None
            moving = self.backend.moving_selection(self.imol, seed, self.options)
            result, hbonds = self.backend.refine(self.imol, moving, stage, self.options)
            self.index += 1  # Acceptance succeeded; never silently rerun an accepted step.
            self.processed_residues.update(seed)
            self.processed_residues.update(_residue_spec(spec) for spec in moving)
            entry = {"window": w, "stage": si, "repeat": repeat, "moving": len(moving),
                     "hbonds": hbonds, "result": result, "cc_before": before, "rotamers": rotamers}
            self.history.append(entry)
            self.message = "%d/%d | %s:%d%s-%d%s | Map #%d, weight %g | %d moving" % (
                self.index, len(self.tasks), seed[0][0], seed[0][1], seed[0][2], seed[-1][1], seed[-1][2],
                stage.map_id, stage.weight, len(moving))
            if rotamers:
                self.message += " | Rotamers %d/%d fitted" % (
                    rotamers["attempted"] - len(rotamers["failed"]), rotamers["attempted"])
                if rotamers["failed"]:
                    print("Multimap rotamer fits unavailable:", rotamers["failed"])
            reasons = []
            if self.options.auto_pause:
                after = self.backend.correlation(self.imol, seed, reference)
                entry["cc_after"] = after
                self.message += " | CC %.3f" % after
                if after < self.options.minimum_mean_cc:
                    reasons.append("mean CC %.3f below %.3f" % (after, self.options.minimum_mean_cc))
                if before - after > self.options.maximum_cc_drop:
                    reasons.append("CC fell %.3f" % (before - after))
                for name, label, value in result[2] or []:
                    if str(name).lower() in ("bonds", "angles", "chirals"):
                        if not math.isfinite(float(value)) or float(value) > self.options.maximum_geometry:
                            reasons.append(str(label))
            if self._requested:
                self._settle(self._requested, self._requested.capitalize() + "; completed edits retained.")
            elif reasons:
                self._settle("paused", self.message + " | CHECK: " + "; ".join(reasons))
            elif self.index == len(self.tasks):
                self._settle("complete", self.message + " | Complete.")
            else:
                self._notify()
        except Exception as error:
            traceback.print_exc()
            self._settle("failed", "Failed: %s. Earlier edits remain; use Coot Undo if needed." % error)
        finally:
            self._busy = False
        return self.state == "running"


class MultimapRefinementRunner:
    """One-shot, low-priority callbacks leave paint time before and after native work."""

    def __init__(self, job, glib=None):
        if glib is None:
            from gi.repository import GLib
            glib = GLib
        self.job, self.glib, self.source = job, glib, None

    def _cancel(self):
        if self.source is not None:
            self.glib.source_remove(self.source)
            self.source = None

    def _queue(self, prepare):
        self.source = self.glib.timeout_add(75, self._dispatch, prepare,
                                           priority=self.glib.PRIORITY_DEFAULT_IDLE)

    def _dispatch(self, prepare):
        self.source = None
        keep = self.job.prepare_step() if prepare else self.job.tick()
        if keep and self.job.state == "running":
            self._queue(not prepare)
        return False

    def start(self):
        self._cancel()
        self.job.start()
        if self.job.state == "running":
            self._queue(True)

    def pause(self):
        self._cancel()
        self.job.pause()

    def stop(self):
        self._cancel()
        self.job.stop()


def start_multimap_refinement(imol, residue_specs, stages, options=None):
    """Start via GLib. Stages are explicit RefinementStage objects, in execution order."""
    job = MultimapRefinementJob(imol, residue_specs, stages, options)
    job.runner = MultimapRefinementRunner(job)
    job.runner.start()
    return job


def multimap_refinement_gui():
    """Map selectors and an ordered stage table; default scope is all polymer chains."""
    from gi.repository import Gtk, GLib, Pango
    backend = CootRefinementBackend()
    api = backend.api
    active = api.active_residue_py()
    if not active:
        api.info_dialog("Choose an active residue first.")
        return
    imol = active[0]
    initial_records = backend.snapshot(imol)
    if not backend.polymer_specs:
        api.info_dialog("The active model has no refinable polymer residues.")
        return
    maps = [i for i in range(api.graphics_n_molecules()) if api.is_valid_map_molecule(i)]
    if not maps:
        api.info_dialog("Load at least one refinement map first.")
        return
    gtk4 = Gtk.get_major_version() >= 4
    window = Gtk.Window(title="Multimap refinement: model #%d" % imol)
    window.set_default_size(650, 480)
    root = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=10)
    for edge in ("top", "bottom", "start", "end"):
        getattr(root, "set_margin_" + edge)(12)
    window.set_child(root) if gtk4 else window.add(root)

    def append(box, child, expand=False):
        box.append(child) if gtk4 else box.pack_start(child, expand, expand, 0)

    def label(text):
        item = Gtk.Label(label=text)
        item.set_xalign(0)
        item.set_wrap(True) if gtk4 else item.set_line_wrap(True)
        return item

    def heading(text):
        item = label("")
        item.set_markup("<b>" + text + "</b>")
        return item

    def child(container, widget):
        container.set_child(widget) if gtk4 else container.add(widget)

    def grid():
        item = Gtk.Grid()
        item.set_row_spacing(6)
        item.set_column_spacing(10)
        return item

    model_label = label("#%d  %s" % (imol, str(api.molecule_name(imol)).split("/")[-1]))
    model_label.set_ellipsize(Pango.EllipsizeMode.MIDDLE)
    append(root, model_label)
    content_scroll = Gtk.ScrolledWindow()
    content_scroll.set_policy(Gtk.PolicyType.NEVER, Gtk.PolicyType.AUTOMATIC)
    content_scroll.set_vexpand(True)
    content_scroll.set_hexpand(True)
    outer = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=12)
    child(content_scroll, outer)
    append(root, content_scroll, True)
    controls = []

    def spin(value, lower, upper, increment=1, digits=0):
        field = Gtk.SpinButton.new_with_range(lower, upper, increment)
        field.set_digits(digits)
        field.set_numeric(True)
        field.set_width_chars(6)
        field.set_value(value)
        return field

    append(outer, heading("Scope"))
    scope_grid = grid()
    scope = Gtk.ComboBoxText()
    polymer_chains = sorted({s[0] for s in backend.polymer_specs})
    chain_counts = chain_residue_counts(initial_records, backend.polymer_specs)
    scope.append("all", "All polymer chains (%d chains; %d polymer, %d non-polymer)" %
                 (len(polymer_chains), sum(chain_counts[name][0] for name in polymer_chains),
                  sum(chain_counts[name][1] for name in polymer_chains)))
    for name in polymer_chains:
        polymer, nonpolymer = chain_counts[name]
        scope.append("chain:" + name, "Chain %s only (%d polymer, %d non-polymer)" %
                     (name or "[blank]", polymer, nonpolymer))
    scope.set_active_id("all")
    scope.set_hexpand(True)
    scope_grid.attach(scope, 0, 0, 4, 1)
    size, step = spin(11, 1, 2000), spin(5, 1, 2000)
    scope_grid.attach(label("Window size"), 0, 1, 1, 1)
    scope_grid.attach(size, 1, 1, 1, 1)
    scope_grid.attach(label("Step size"), 2, 1, 1, 1)
    scope_grid.attach(step, 3, 1, 1, 1)
    append(outer, scope_grid)
    controls.append(scope_grid)

    append(outer, heading("Refinement stages"))
    stage_grid = grid()
    append(outer, stage_grid)
    controls.append(stage_grid)
    stage_rows = []
    selected_map = api.imol_refinement_map()
    if selected_map not in maps:
        selected_map = maps[0]

    def populate_maps(combo, chosen):
        combo.remove_all()
        for map_id in maps:
            name = str(api.molecule_name(map_id)).split("/")[-1]
            combo.append(str(map_id), "#%d  %s" % (map_id, name))
        combo.set_active_id(str(chosen))
        if combo.get_active_id() is not None:
            combo.set_tooltip_text(str(api.molecule_name(int(combo.get_active_id()))))
        else:
            combo.set_tooltip_text("Choose a map; the previous map may have been closed.")

    def redraw_stages():
        if gtk4:
            while stage_grid.get_first_child():
                stage_grid.remove(stage_grid.get_first_child())
        else:
            for widget in stage_grid.get_children():
                stage_grid.remove(widget)
        for column, text in enumerate(("", "Map", "Weight", "Repeats", "Autofit\nrotamers", "")):
            stage_grid.attach(label(text), column, 0, 1, 1)
        for index, row in enumerate(stage_rows):
            stage_grid.attach(label(str(index+1)), 0, index+1, 1, 1)
            for column, key in enumerate(("map", "weight", "repeats", "autofit_rotamers"), 1):
                stage_grid.attach(row[key], column, index+1, 1, 1)
            actions = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=2)
            for text, change, tooltip in (("^", -1, "Move stage earlier"),
                                          ("v", 1, "Move stage later"),
                                          ("-", 0, "Remove stage")):
                button = Gtk.Button(label=text)
                button.set_tooltip_text(tooltip)
                button.set_sensitive((len(stage_rows) > 1 if not change else
                                      0 <= index + change < len(stage_rows)))
                def edit(_button, item=row, delta=change):
                    current = stage_rows.index(item)
                    if delta:
                        stage_rows[current], stage_rows[current+delta] = stage_rows[current+delta], stage_rows[current]
                    else:
                        stage_rows.remove(item)
                    redraw_stages()
                button.connect("clicked", edit)
                append(actions, button)
            stage_grid.attach(actions, 5, index+1, 1, 1)
        if not gtk4:
            stage_grid.show_all()

    def add_stage(_button=None):
        chosen = stage_rows[-1]["map"].get_active_id() if stage_rows else str(selected_map)
        weight = stage_rows[-1]["weight"].get_value() if stage_rows else api.matrix_state()
        combo = Gtk.ComboBoxText()
        combo.set_hexpand(True)
        for renderer in combo.get_cells():
            renderer.set_property("ellipsize", Pango.EllipsizeMode.MIDDLE)
            renderer.set_property("max-width-chars", 32)
        populate_maps(combo, chosen)
        def map_changed(item):
            value = item.get_active_id()
            item.set_tooltip_text(str(api.molecule_name(int(value))) if value else "Choose a map")
        combo.connect("changed", map_changed)
        autofit = Gtk.CheckButton()
        autofit.set_tooltip_text("Autofit this range's protein sidechains against this map before each refinement. Uses current rotamer/backrub mode.")
        autofit.set_halign(Gtk.Align.CENTER)
        stage_rows.append({"map": combo, "weight": spin(weight, 0.01, 100000, 1, 2),
                           "repeats": spin(1, 1, 100), "autofit_rotamers": autofit})
        redraw_stages()

    def refresh_maps(_button=None):
        maps[:] = [i for i in range(api.graphics_n_molecules()) if api.is_valid_map_molecule(i)]
        for row in stage_rows:
            populate_maps(row["map"], row["map"].get_active_id())

    add_stage()
    stage_toolbar = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=6)
    for title, callback in (("Add stage", add_stage), ("Refresh maps", refresh_maps)):
        item = Gtk.Button(label=title)
        item.connect("clicked", callback)
        append(stage_toolbar, item)
    append(outer, stage_toolbar)
    controls.append(stage_toolbar)
    hb = Gtk.CheckButton(label="Temporary backbone H-bond restraints")
    guard = Gtk.CheckButton(label="Auto-pause on poor fit or geometry")
    for item in (hb, guard):
        append(outer, item)
        controls.append(item)
    advanced = Gtk.Expander(label="Options")
    advanced_box = grid()
    child(advanced, advanced_box)
    radius, padding, limit = spin(4, 0, 30, 0.5, 1), spin(1, 0, 10), spin(0, 0, 10000)
    limit.set_tooltip_text("0 = automatic: at least 500 residues or 5 times the window size. A positive value sets an explicit cap.")
    min_cc, cc_drop, geometry = spin(0.2, -1, 1, 0.05, 2), spin(0.1, 0, 2, 0.05, 2), spin(5, 0.1, 100, 0.5, 1)
    for index, (text, field) in enumerate((("Contact shell (A)", radius), ("Neighbour padding", padding),
                                          ("Moving limit (0 = auto)", limit), ("Minimum mean CC", min_cc),
                                          ("Maximum CC drop", cc_drop), ("Maximum geometry distortion", geometry))):
        advanced_box.attach(label(text), 0, index, 1, 1)
        advanced_box.attach(field, 1, index, 1, 1)
    controls.append(advanced_box)
    append(outer, advanced)
    append(outer, label("Each range: all stages top to bottom, then the next range.\nCC uses the last map. Pause/Stop act between steps; completed edits remain."))
    changes = Gtk.Expander(label="Changes from starting model")
    changes.set_sensitive(False)
    changes_box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6)
    child(changes, changes_box)
    changes_row = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=6)
    change_metric = Gtk.ComboBoxText()
    for key, text in (("ca", "C-alpha RMSD (sliding window)"), ("chi", "Rotamer torsion changes"),
                      ("backbone", "Backbone dihedral changes")):
        change_metric.append(key, text)
    change_metric.set_active_id("ca")
    change_metric.set_hexpand(True)
    change_window = spin(5, 1, 101, 2)
    change_window.set_tooltip_text("Odd-sized, centred C-alpha window; truncated at chain breaks and termini.")
    show_changes = Gtk.Button(label="Colour model")
    for item in (change_metric, label("Window"), change_window, show_changes):
        append(changes_row, item)
    append(changes_box, changes_row)
    change_description = label("")
    change_legend = label("")
    change_legend.set_selectable(True)
    change_summary = label("")
    for item in (change_description, change_legend, change_summary):
        append(changes_box, item)

    def change_selection(_combo=None):
        metric = change_metric.get_active_id()
        change_window.set_sensitive(metric == "ca")
        change_description.set_text({"ca": "Map-frame displacement; no superposition.",
                                     "chi": "Largest wrapped chi change; symmetric terminal torsions folded.",
                                     "backbone": "Largest wrapped phi, psi or omega change."}[metric])
        change_legend.set_text(change_colour_legend(metric))

    change_metric.connect("changed", change_selection)
    change_selection()
    append(outer, changes)
    progress = Gtk.ProgressBar()
    progress.set_show_text(True)
    append(root, progress)
    status = label("Ready.")
    append(root, status)
    buttons = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=6)
    append(root, buttons)
    state = {"job": None, "runner": None, "comparison_shown": False}

    def visualize(_button):
        try:
            change_window.update()
            metric = change_metric.get_active_id()
            values = state["job"].visualize_changes(metric, change_window.get_value_as_int() if metric == "ca" else 5)
            scorable = [v for v in values.values() if v is not None]
            if scorable:
                text = "%d/%d residues compared | Mean %.2f, max %.2f %s" % (
                    len(scorable), len(values), sum(scorable)/len(scorable), max(scorable), _CHANGE_METRICS[metric][1])
            else:
                text = "0/%d residues compared | No processed residues with comparable data." % len(values)
            change_summary.set_text(text)
            api.add_status_bar_text(text)
        except Exception as error:
            change_summary.set_text(str(error))
            traceback.print_exc()

    show_changes.connect("clicked", visualize)

    def update(job):
        progress.set_fraction(job.index / len(job.tasks))
        progress.set_text("%d / %d steps complete" % (job.index, len(job.tasks)))
        status.set_text(job.message)
        busy = job.state in ("running", "paused")
        for control in controls:
            control.set_sensitive(not busy)
        run_button.set_sensitive(job.state != "running")
        run_button.set_label("Resume" if job.state == "paused" else "Start")
        pause_button.set_sensitive(job.state == "running")
        stop_button.set_sensitive(busy)
        available = job.state in ("paused", "complete", "stopped", "failed") and job.initial_coordinates is not None
        changes.set_sensitive(available)
        if available and not state["comparison_shown"]:
            changes.set_expanded(True)
            state["comparison_shown"] = True
        api.add_status_bar_text(job.message)

    def run(_button):
        try:
            job = state["job"]
            if job is None or job.state != "paused":
                backend.snapshot(imol)
                chosen_scope = scope.get_active_id()
                specs = sorted(s for s in backend.polymer_specs if chosen_scope == "all" or
                               (chosen_scope and s[0] == chosen_scope[6:]))
                stages = []
                for index, row in enumerate(stage_rows):
                    map_id = row["map"].get_active_id()
                    if map_id is None or not api.is_valid_map_molecule(int(map_id)):
                        raise ValueError("Choose an available map for stage %d." % (index+1))
                    row["weight"].update()
                    row["repeats"].update()
                    stages.append(RefinementStage(int(map_id), row["weight"].get_value(),
                                                  row["repeats"].get_value_as_int(),
                                                  autofit_rotamers=row["autofit_rotamers"].get_active()))
                for field in (size, step, radius, padding, limit, min_cc, cc_drop, geometry):
                    field.update()
                contact = radius.get_value()
                options = RefinementOptions(window_size=size.get_value_as_int(), step_size=step.get_value_as_int(),
                                            contact_radius=contact, neighbour_padding=padding.get_value_as_int() if contact else 0,
                                            max_moving_residues=limit.get_value_as_int(),
                                            hydrogen_bonds=hb.get_active(), auto_pause=guard.get_active(),
                                            minimum_mean_cc=min_cc.get_value(), maximum_cc_drop=cc_drop.get_value(),
                                            maximum_geometry=geometry.get_value(),
                                            order="window")
                refinement_windows(specs, options.window_size, options.step_size)
                backend._check_idle()
                if not api.is_valid_map_molecule(api.imol_refinement_map()):
                    api.set_imol_refinement_map(stages[0].map_id)
                job = MultimapRefinementJob(imol, specs, stages, options, backend, update)
                state["job"] = job
                state["comparison_shown"] = False
                change_summary.set_text("")
                state["runner"] = MultimapRefinementRunner(job, GLib)
            state["runner"].start()
        except Exception as error:
            status.set_text(str(error))
            traceback.print_exc()

    def pause(_button):
        if state["runner"]:
            state["runner"].pause()

    def stop(_button=None):
        if state["runner"]:
            state["runner"].stop()

    def button(text, callback):
        item = Gtk.Button(label=text)
        item.connect("clicked", callback)
        append(buttons, item)
        return item

    run_button = button("Start", run)
    pause_button = button("Pause", pause)
    stop_button = button("Stop", stop)
    pause_button.set_sensitive(False)
    stop_button.set_sensitive(False)
    # Retain references for native widget smoke tests and interactive inspection.
    window.multimap_widgets = {"scope": scope, "stages": stage_rows, "add_stage": add_stage,
                               "refresh_maps": refresh_maps, "start": run_button, "status": status,
                               "progress": progress, "pause": pause_button,
                               "stop": stop_button, "changes": changes, "change_metric": change_metric,
                               "change_window": change_window, "show_changes": show_changes,
                               "change_summary": change_summary}
    window.connect("close-request" if gtk4 else "delete-event", lambda *_: (stop(), False)[1])
    window.present() if gtk4 else window.show_all()
    return window
'''


def _make_multimap_user_defined_colours(module):
  colours = _make_default_user_defined_colours()
  return colours + module._standalone_change_colour_table()[len(colours):]


def _multimap_refinement_module():
  cached = getattr(coot, "_trimmings_multimap_refinement", None)
  if cached is not None:
    return cached
  from types import ModuleType
  name = "_coot_trimmings_multimap_refinement"
  module = ModuleType(name)
  # Dataclasses resolve type metadata through the registered module.
  sys.modules[name] = module
  try:
    exec(compile(_MULTIMAP_REFINEMENT_SOURCE, "<coot-trimmings-multimap>", "exec"),
         module.__dict__)
    module._standalone_change_colour_table = module._change_colour_table
    module._change_colour_table = lambda: _make_multimap_user_defined_colours(module)

    class TrimmingsMultimapBackend(module.CootRefinementBackend):
      def validate(self, *args, **kwargs):
        if globals().get("_USER_DEFINED_COLOUR_JOB_ACTIVE"):
          raise RuntimeError("Wait for model colouring to finish before starting refinement.")
        return super().validate(*args, **kwargs)

      def autofit_rotamers(self, imol, seed, stage):
        try:
          return super().autofit_rotamers(imol, seed, stage)
        finally:
          _mark_model_changed(imol)

      def refine(self, imol, specs, stage, options):
        try:
          return super().refine(imol, specs, stage, options)
        finally:
          # Failed native operations can also leave partial coordinate edits.
          _mark_model_changed(imol)

      def colour_changes(self, imol, records, values, metric):
        result = super().colour_changes(imol, records, values, metric)
        representation = (graphics_to_user_defined_atom_colours_all_atoms_representation
                          if metric == "chi" else
                          graphics_to_user_defined_atom_colours_representation)
        sequence = _representation_sequence()
        if representation in sequence:
          cycle_rep_flag[imol] = sequence.index(representation)
        return result

    module.CootRefinementBackend = TrimmingsMultimapBackend
  except Exception:
    if sys.modules.get(name) is module:
      del sys.modules[name]
    raise
  # Keep paused/running jobs intact even if the startup script is reloaded.
  coot._trimmings_multimap_refinement = module
  return module


def multimap_refinement_gui():
  """Open stepped multimap refinement for the active model."""
  try:
    return _multimap_refinement_module().multimap_refinement_gui()
  except Exception as error:
    traceback.print_exc()
    info_dialog("Multimap stepped refinement: " + str(error))
    return None



#Refine triple (Paul)
def key_binding_refine_triple():
    active_atom = active_residue()
    if not active_atom:
       add_status_bar_text("No active atom")
       return None
    else:
       imol = active_atom[0]
       residue_spec = _coot_residue_spec_from_spec(atom_spec_to_residue_spec(active_atom))
       if residue_spec is None:
         add_status_bar_text("Could not resolve the active residue for triple refinement")
         return None
       N_terminal_residue_spec = [residue_spec[0], residue_spec[1]-1, residue_spec[2]]
       C_terminal_residue_spec = [residue_spec[0], residue_spec[1]+1, residue_spec[2]]
       spec_list = [N_terminal_residue_spec, list(residue_spec), C_terminal_residue_spec]
       refine_residues(imol, spec_list)
    
#Cycle symmetry represntation mode forward/back
cycle_symm_flag={0:0}
def cycle_symm_up(mol_id,flag):
  global cycle_symm_flag
  cycle_symm_flag[mol_id]=flag
  if cycle_symm_flag[mol_id]==0:
    get_show_symmetry()
    set_show_symmetry_master(1)
    symmetry_as_calphas(mol_id,1)
    cycle_symm_flag[mol_id]=1
  elif cycle_symm_flag[mol_id]==1:
    get_show_symmetry()
    set_show_symmetry_master(1)
    symmetry_as_calphas(mol_id,0)
    set_symmetry_whole_chain(mol_id,1)
    cycle_symm_flag[mol_id]=2
  elif cycle_symm_flag[mol_id]==2:
    get_show_symmetry()
    set_show_symmetry_master(1)
    symmetry_as_calphas(mol_id,0)
    set_symmetry_whole_chain(mol_id,0)
    cycle_symm_flag[mol_id]=3
  elif cycle_symm_flag[mol_id]==3:
    get_show_symmetry()
    set_show_symmetry_master(0)
    cycle_symm_flag[mol_id]=0
def cycle_symm_down(mol_id,flag):
  global cycle_symm_flag
  cycle_symm_flag[mol_id]=flag
  if cycle_symm_flag[mol_id]==3:
    get_show_symmetry()
    set_show_symmetry_master(1)
    symmetry_as_calphas(mol_id,0)
    set_symmetry_whole_chain(mol_id,1)
    cycle_symm_flag[mol_id]=2
  elif cycle_symm_flag[mol_id]==2:
    get_show_symmetry()
    set_show_symmetry_master(1)
    symmetry_as_calphas(mol_id,1)
    cycle_symm_flag[mol_id]=1
  elif cycle_symm_flag[mol_id]==1:
    get_show_symmetry()
    set_show_symmetry_master(0)
    cycle_symm_flag[mol_id]=0
  elif cycle_symm_flag[mol_id]==0:
    get_show_symmetry()
    set_show_symmetry_master(1)
    symmetry_as_calphas(mol_id,0)
    set_symmetry_whole_chain(mol_id,0)
    cycle_symm_flag[mol_id]=3


def cycle_symm_up_current():
  residue = _active_residue_or_status()
  if not residue:
    return None
  mol_id = residue[0]
  cycle_symm_up(mol_id, cycle_symm_flag.get(mol_id, 0))


def cycle_symm_down_current():
  residue = _active_residue_or_status()
  if not residue:
    return None
  mol_id = residue[0]
  cycle_symm_down(mol_id, cycle_symm_flag.get(mol_id, 0))
    
def undo_visible():
  residue = _active_residue_or_status()
  if not residue:
    return None
  set_undo_molecule(residue[0])
  apply_undo()
  _mark_model_changed(residue[0])
  
def redo_visible():
  residue = _active_residue_or_status()
  if not residue:
    return None
  set_undo_molecule(residue[0])
  apply_redo()
  _mark_model_changed(residue[0])

#Cycle rotamers for active residue
rotamer_number=0
def cycle_rotamers():
  global rotamer_number
  residue = _active_residue_or_status()
  if not residue:
    return None
  mol_id=residue[0]
  ch_id=residue[1]
  res_here=residue[2]
  ins_code=""
  alt_conf=""
  n_rots=n_rotamers(mol_id, ch_id, res_here, ins_code)-1
  with _grouped_model_edit(mol_id):
    update_go_to_atom_from_current_position()
    if rotamer_number>=n_rots:
      rotamer_number=0
      set_residue_to_rotamer_number(mol_id,ch_id,res_here,ins_code,alt_conf,rotamer_number)
    else:
      rotamer_number=rotamer_number+1
      set_residue_to_rotamer_number(mol_id,ch_id,res_here,ins_code,alt_conf,rotamer_number)
  
#Toggle environment distances
def toggle_env_dist():
  if show_environment_distances_state()==1:
    set_show_environment_distances(0)
  else:
    set_show_environment_distances(1)
    
#Check if residue is in a polymer by comparing resname to list
_SEGMENT_POLYMER_NAMES = set("A C T G U DA DC DT DG DU ALA UNK ARG ASN ASP CYS GLU GLN GLY HIS ILE LEU LYS MET MSE PHE PRO SER THR TRP TYR VAL".split())


def _segment_chain_snapshot(mol_id, ch_id, model_snapshot=None):
  if model_snapshot is None and callable(getattr(coot, "python_representation_kk", None)):
    model_snapshot = _model_coordinate_snapshot(mol_id)
  rows = []
  if model_snapshot is not None:
    entries = [(number, ins, record["name"], record["atoms"])
               for (chain, number, ins), record in model_snapshot.items() if chain == ch_id]
  else:
    entries = []
    for sn in range(chain_n_residues(ch_id, mol_id)):
      number = seqnum_from_serial_number(mol_id, ch_id, sn)
      ins = insertion_code_from_serial_number(mol_id, ch_id, sn)
      entries.append((number, ins, resname_from_serial_number(mol_id, ch_id, sn),
                      residue_info_py(mol_id, ch_id, number, ins) or []))
  for number, ins, name, atoms in entries:
    xyz = {}
    for atom in sorted(atoms, key=lambda a: (bool(a[0][1]), -float(a[1][0]))):
      xyz.setdefault(atom[0][0].strip(), atom[2])
    rows.append({"number": number, "ins": ins,
                 "name": name, "xyz": xyz})
  return rows


def _segment_rows_connected(left, right):
  if left["name"] not in _SEGMENT_POLYMER_NAMES or right["name"] not in _SEGMENT_POLYMER_NAMES:
    return False
  for atom1, atom2 in (("C", "N"), ("O3'", "P"), ("O3*", "P")):
    a, b = left["xyz"].get(atom1), right["xyz"].get(atom2)
    if a is not None and b is not None:
      return 0.25 < sum((x-y)**2 for x, y in zip(a, b)) <= 2.2**2
  return False


def _segment_end_for_residue(mol_id, ch_id, resno, first):
  rows = _segment_chain_snapshot(mol_id, ch_id)
  matches = [i for i, row in enumerate(rows) if row["number"] == resno]
  if len(matches) != 1 or rows[matches[0]]["name"] not in _SEGMENT_POLYMER_NAMES:
    return -1
  index = matches[0]
  if first:
    while index > 0 and _segment_rows_connected(rows[index-1], rows[index]):
      index -= 1
  else:
    while index+1 < len(rows) and _segment_rows_connected(rows[index], rows[index+1]):
      index += 1
  return rows[index]["number"]


def is_polymer_residue(mol_id, ch_id, sn):
  if not 0 <= sn < chain_n_residues(ch_id, mol_id):
    return 0
  return int(resname_from_serial_number(mol_id, ch_id, sn) in _SEGMENT_POLYMER_NAMES)
    
#Return last residue in polymer
def last_polymer_residue(mol_id, ch_id):
  for sn in reversed(range(chain_n_residues(ch_id, mol_id))):
    if is_polymer_residue(mol_id, ch_id, sn):
      return seqnum_from_serial_number(mol_id, ch_id, sn)
  return -1
    
def first_polymer_residue(mol_id, ch_id):
  for sn in range(chain_n_residues(ch_id, mol_id)):
    if is_polymer_residue(mol_id, ch_id, sn):
      return seqnum_from_serial_number(mol_id, ch_id, sn)
  return -1


#Return last res (polymer or not)
def last_residue(mol_id,ch_id):
          n=chain_n_residues(ch_id,mol_id)-1
          result=seqnum_from_serial_number(mol_id,"%s"%(ch_id),n)
          return result

#Return first res (polymer or not)
def first_residue(mol_id,ch_id):
          result=seqnum_from_serial_number(mol_id,"%s"%(ch_id),0)
          return result
          
#Check if chain is polymer
def is_polymer(mol_id,ch_id):
        a=is_protein_chain_p(mol_id,"%s"%(ch_id))
        b=is_nucleotide_chain_p(mol_id,"%s"%(ch_id))
        if (a==1) or (b==1):
                result=1
        else:
                result=0
        return result
    
#Check if residue is last in polymer
def is_last_polymer_residue_sn(mol_id, ch_id, sn):
  return int(is_polymer_residue(mol_id, ch_id, sn) and not is_polymer_residue(mol_id, ch_id, sn+1))

#get serial number from resnum
def get_sn_from_resno(mol_id, ch_id, resno):
  matches = [sn for sn in range(chain_n_residues(ch_id, mol_id))
             if seqnum_from_serial_number(mol_id, ch_id, sn) == resno]
  return matches[0] if len(matches) == 1 else -1

def get_sn_from_resno_alt(mol_id,ch_id,resno):
  sn=0
  sn_max=chain_n_residues(ch_id,mol_id)-1
  sn_dict={}
  while sn<=sn_max:
    resno_here=seqnum_from_serial_number(mol_id,ch_id,sn)
    sn_dict[resno_here]=sn
    sn=sn+1
  try:
    sn_out=sn_dict[resno]
    return sn_out
  except KeyError:
    return -1

#check if res is at C-term side of break in mid chain
def is_term_type_mc(mol_id, ch_id, resno):
  return is_term_type_mc_sn(mol_id, ch_id, get_sn_from_resno(mol_id, ch_id, resno))
    

#check if res is at C-term side of break in mid chain
def is_term_type_mc_sn(mol_id, ch_id, sn):
  rows = _segment_chain_snapshot(mol_id, ch_id)
  return int(0 <= sn < len(rows)-1 and rows[sn]["name"] in _SEGMENT_POLYMER_NAMES
             and rows[sn+1]["name"] in _SEGMENT_POLYMER_NAMES
             and not _segment_rows_connected(rows[sn], rows[sn+1]))
    
#check if res is at N-term side of break in mid chain
def is_term_type_mn(mol_id, ch_id, resno):
  return is_term_type_mn_sn(mol_id, ch_id, get_sn_from_resno(mol_id, ch_id, resno))
    
def is_term_type_mn_sn(mol_id, ch_id, sn):
  return is_term_type_mc_sn(mol_id, ch_id, sn-1)
  
        
def first_residue_in_seg(mol_id, ch_id, resno):
  return _segment_end_for_residue(mol_id, ch_id, resno, first=True)
  
def last_residue_in_seg(mol_id, ch_id, resno):
  return _segment_end_for_residue(mol_id, ch_id, resno, first=False)


#Get serial number of active residue
def sn_of_active_res():
  chain_context = _active_chain_context_or_status()
  if not chain_context:
    return None
  mol_id = chain_context["mol_id"]
  ch_id = chain_context["chain_id"]
  resn_to_match = chain_context["resno"]
  for sn in _iter_chain_serial_numbers(mol_id, ch_id):
    if seqnum_from_serial_number(mol_id, "%s" % ch_id, sn) == resn_to_match:
      return sn
  return None
    
#Get monomer and delete hydrogens
def get_monomer_no_H(mon):
  existing = set(model_molecule_list())
  mol_id = get_monomer(mon)
  if mol_id in existing or not valid_model_molecule_qm(mol_id):
    info_dialog(f"Could not create monomer {mon}; existing models were not changed.")
    return -1
  delete_hydrogens(mol_id)
  return mol_id
    
#Return list of segments in active mol
def segment_list(mol_id):
  return [[mol_id, specs[0][0], specs[0][1], specs[-1][1]]
          for specs in _model_segment_specs(mol_id)]


def _model_segment_specs(mol_id):
  snapshot = _model_coordinate_snapshot(mol_id)
  result = []
  for chain in dict.fromkeys(spec[0] for spec in snapshot):
    previous = None
    for row in _segment_chain_snapshot(mol_id, chain, model_snapshot=snapshot):
      if row["name"] not in _SEGMENT_POLYMER_NAMES:
        previous = None
        continue
      if previous is None or not _segment_rows_connected(previous, row):
        result.append([])
      result[-1].append((chain, row["number"], row["ins"]))
      previous = row
  return result
  
def segment_list_chain(mol_id, ch_id):
  rows = _segment_chain_snapshot(mol_id, ch_id)
  result = []
  previous = None
  for row in rows:
    if row["name"] not in _SEGMENT_POLYMER_NAMES:
      previous = None
      continue
    if previous is None or not _segment_rows_connected(previous, row):
      result.append([mol_id, ch_id, row["number"], row["number"]])
    else:
      result[-1][3] = row["number"]
    previous = row
  return result
  
    
#Local cylinder refinement around the active residue:
# main polymer range +/- 5 residues, contact-expanded secondary ranges,
# then local polymer gap-filling/pruning while keeping non-polymer contacts.
def auto_refine(primary_specs=None, mol_id=None, secondary_structure=None):
  if imol_refinement_map()==-1:
    info_dialog("You must set a refinement map!")
    return None
  if coot.refinement_already_ongoing_p():
    info_dialog("Accept or reject the current refinement first.")
    return None
  if primary_specs is None:
    residue = _active_residue_or_status()
    if not residue:
      return None
    mol_id=residue[0]
  else:
    primary_specs = [tuple(spec) for spec in primary_specs]
    if not primary_specs or mol_id not in model_molecule_list():
      raise ValueError("Select a range in an available model.")
    residue = [mol_id, *primary_specs[0]]
  ch_id=residue[1]
  resno=residue[2]
  half_window = 5
  contact_radius = 4.0
  max_gap_to_fill = 4

  def nearest_polymer_residue(molecule_id, chain_id, residue_no, ins_code):
    normalized_active_spec = _coot_residue_spec_from_spec([chain_id, residue_no, ins_code])
    if normalized_active_spec is None:
      return None
    active_spec = list(normalized_active_spec)
    if is_polymer_residue_spec(molecule_id, normalized_active_spec):
      return active_spec

    for search_radius in (4.0, 6.0, 8.0, 10.0):
      nearby_specs = residues_near_residue(molecule_id, active_spec, search_radius)
      polymer_specs = []
      for spec in nearby_specs or []:
        normalized_spec = _coot_residue_spec_from_spec(spec)
        if normalized_spec is None:
          continue
        if is_polymer_residue_spec(molecule_id, normalized_spec):
          polymer_specs.append(list(normalized_spec))
      if polymer_specs:
        return polymer_specs[0]

    best_residue = None
    best_distance = None
    for _mol_id, _chain_id, seg_start, seg_end in segment_list_chain(molecule_id, chain_id):
      if residue_no < seg_start:
        candidate_residues = [seg_start]
      elif residue_no > seg_end:
        candidate_residues = [seg_end]
      else:
        candidate_residues = [residue_no]
      for candidate in candidate_residues:
        distance = abs(candidate - residue_no)
        if best_distance is None or distance < best_distance:
          best_distance = distance
          best_residue = candidate
    return best_residue

  def is_polymer_residue_spec(molecule_id, residue_spec):
    chain_id, residue_no, ins_code = residue_spec
    return _residue_is_polymer(molecule_id, chain_id, residue_no, ins_code)

  def fill_short_polymer_gaps(molecule_id, polymer_specs):
    expanded_specs = set(polymer_specs)
    polymer_by_chain = {}
    for chain_id, residue_no, ins_code in polymer_specs:
      if ins_code == "":
        polymer_by_chain.setdefault(chain_id, set()).add(residue_no)
    for chain_id, residue_numbers in polymer_by_chain.items():
      for _mol_id, segment_chain_id, segment_start, segment_end in segment_list_chain(molecule_id, chain_id):
        residues_in_segment = sorted(
          residue_no for residue_no in residue_numbers if segment_start <= residue_no <= segment_end
        )
        for index in range(len(residues_in_segment)-1):
          left_residue = residues_in_segment[index]
          right_residue = residues_in_segment[index+1]
          gap_size = right_residue - left_residue - 1
          if gap_size <= 0 or gap_size > max_gap_to_fill:
            continue
          for fill_residue in range(left_residue+1, right_residue):
            if residue_exists_qm(molecule_id, segment_chain_id, fill_residue, ""):
              expanded_specs.add((segment_chain_id, fill_residue, ""))
    return expanded_specs

  def expand_polymer_contact_windows(molecule_id, polymer_specs, window_size):
    expanded_specs = set()
    polymer_by_chain = {}
    for chain_id, residue_no, ins_code in polymer_specs:
      if ins_code == "":
        polymer_by_chain.setdefault(chain_id, set()).add(residue_no)
      else:
        expanded_specs.add((chain_id, residue_no, ins_code))
    for chain_id, residue_numbers in polymer_by_chain.items():
      fpr = first_polymer_residue(molecule_id, chain_id)
      lpr = last_polymer_residue(molecule_id, chain_id)
      if fpr == -10000 or lpr == -10000:
        continue
      for residue_no in residue_numbers:
        res_start = max(fpr, residue_no - window_size)
        res_end = min(lpr, residue_no + window_size)
        for resn in range(res_start, res_end + 1):
          if residue_exists_qm(molecule_id, chain_id, resn, ""):
            expanded_specs.add((chain_id, resn, ""))
    return expanded_specs

  def prune_small_polymer_fragments(molecule_id, polymer_specs, minimum_fragment_size):
    kept_specs = set()
    polymer_by_chain = {}
    for chain_id, residue_no, ins_code in polymer_specs:
      if ins_code == "":
        polymer_by_chain.setdefault(chain_id, set()).add(residue_no)
    for chain_id, residue_numbers in polymer_by_chain.items():
      for _mol_id, segment_chain_id, segment_start, segment_end in segment_list_chain(molecule_id, chain_id):
        residues_in_segment = sorted(
          residue_no for residue_no in residue_numbers if segment_start <= residue_no <= segment_end
        )
        if not residues_in_segment:
          continue
        cluster = [residues_in_segment[0]]
        for residue_no in residues_in_segment[1:]:
          if residue_no == cluster[-1] + 1:
            cluster.append(residue_no)
          else:
            if len(cluster) >= minimum_fragment_size:
              for cluster_residue in cluster:
                kept_specs.add((segment_chain_id, cluster_residue, ""))
            cluster = [residue_no]
        if len(cluster) >= minimum_fragment_size:
          for cluster_residue in cluster:
            kept_specs.add((segment_chain_id, cluster_residue, ""))
    return kept_specs

  active_ins_code = residue[3]
  main_range_specs = set()
  seed_residue_specs = []

  if primary_specs is not None:
    centre_residue_spec = primary_specs[0]
    if any(not residue_exists_qm(mol_id, *spec) or not _residue_is_polymer(mol_id, *spec)
           for spec in primary_specs):
      raise ValueError("The primary range must contain existing polymer residues only.")
    main_range_specs.update(primary_specs)
    seed_residue_specs.extend(list(spec) for spec in primary_specs)
  elif _residue_is_polymer(mol_id, ch_id, resno, active_ins_code):
    centre_residue_spec = (ch_id, resno, active_ins_code)
    fpr=first_polymer_residue(mol_id,ch_id)
    lpr=last_polymer_residue(mol_id,ch_id)
    res_start=max(fpr, resno-half_window)
    res_end=min(lpr, resno+half_window)

    for resn in range(res_start, res_end+1):
      if residue_exists_qm(mol_id,ch_id,resn,""):
        spec=(ch_id,resn,"")
        seed_residue_specs.append([ch_id,resn,""])
        main_range_specs.add(spec)
  else:
    centre_residue_spec = (ch_id, resno, active_ins_code)
    if not residue_exists_qm(mol_id, ch_id, resno, active_ins_code):
      add_status_bar_text("No active residue for local cylinder refinement")
      return None
    seed_residue_specs.append([ch_id, resno, active_ins_code])
    main_range_specs.add(centre_residue_spec)

  secondary_residue_specs = set()
  direct_non_polymer_contact_specs = set()
  for spec in seed_residue_specs:
    for nearby_residue in residues_near_residue(mol_id, spec, contact_radius):
      nearby_spec = _coot_residue_spec_from_spec(nearby_residue)
      if nearby_spec is not None:
        if nearby_spec not in main_range_specs:
          secondary_residue_specs.add(nearby_spec)
          if not is_polymer_residue_spec(mol_id, nearby_spec):
            direct_non_polymer_contact_specs.add(nearby_spec)

  secondary_polymer_specs = {spec for spec in secondary_residue_specs if is_polymer_residue_spec(mol_id, spec)}
  secondary_non_polymer_specs = (secondary_residue_specs - secondary_polymer_specs) | direct_non_polymer_contact_specs
  if not _residue_is_polymer(mol_id, ch_id, resno, active_ins_code):
    secondary_polymer_specs = expand_polymer_contact_windows(mol_id, secondary_polymer_specs, 1)
  secondary_polymer_specs = fill_short_polymer_gaps(mol_id, secondary_polymer_specs)
  secondary_polymer_specs = prune_small_polymer_fragments(mol_id, secondary_polymer_specs, 3)

  selected_residue_specs = main_range_specs | secondary_polymer_specs | secondary_non_polymer_specs

  residue_list = [list(spec) for spec in sorted(selected_residue_specs, key=lambda spec: (spec[0], spec[1], spec[2]))]
  return _start_cylinder_refinement(mol_id, residue_list, sorted(main_range_specs), secondary_structure)


def _cylinder_secondary_restraints(snapshot, primary, kind):
  """Coot's local helix/strand targets, confined to connected primary residues."""
  if kind not in ("alpha helix", "beta strand"):
    raise ValueError("Unknown cylinder secondary-structure mode.")
  if len(primary) < (4 if kind == "alpha helix" else 3):
    raise ValueError("Select at least four residues for a helix, or three for a strand.")
  atoms = {}
  for spec in primary:
    groups = {}
    for atom in snapshot[spec]["atoms"]:
      if float(atom[1][0]) > 0:
        groups.setdefault(atom[0][0].strip(), []).append(atom)
    if not {"N", "CA", "C", "O"} <= groups.keys():
      raise ValueError("Secondary-structure refinement needs complete protein N/CA/C/O atoms in the primary range.")
    atoms[spec] = groups
  def compatible(specs, names):
    for selected in product(*(atoms[spec][name] for spec, name in zip(specs, names))):
      if len({a[0][1] for a in selected if a[0][1]}) <= 1:
        yield [[*spec, atom[0][0], atom[0][1]] for spec, atom in zip(specs, selected)]
  runs = []
  for spec in primary:
    connected = False
    if runs:
      previous = runs[-1][-1]
      if spec[0] == previous[0] and spec[1] - previous[1] in (0, 1):
        connected = any(0.25 < _distance_sq(a[2], b[2]) <= 2.2**2
                        for a in atoms[previous]["C"] for b in atoms[spec]["N"]
                        if not a[0][1] or not b[0][1] or a[0][1] == b[0][1])
    if not connected:
      runs.append([])
    runs[-1].append(spec)
  restraints = []
  for run in runs:
    for i, spec in enumerate(run):
      if kind == "alpha helix":
        for step, on_target, oo_target in ((3, 3.181, 4.92), (4, 2.919, 6.16)):
          if i + step >= len(run):
            continue
          other = run[i + step]
          if snapshot[other]["name"] != "PRO":
            restraints.extend(("bond", pair, on_target, 0.04)
                              for pair in compatible((spec, other), ("O", "N")))
          restraints.extend(("bond", pair, oo_target, 0.07)
                            for pair in compatible((spec, other), ("O", "O")))
      else:
        if i + 1 < len(run):
          restraints.extend(("bond", pair, 4.64, 0.08)
                            for pair in compatible(run[i:i+2], ("O", "O")))
        if i + 2 < len(run):
          for name, target in (("O", 98.0), ("CA", 120.0)):
            restraints.extend(("angle", triple, target, 0.5)
                              for triple in compatible(run[i:i+3], (name, name, name)))
  if not restraints:
    raise ValueError("No connected primary residues support these secondary-structure restraints.")
  return restraints


def _cylinder_restraint_key(kind, specs):
  atoms = tuple((str(s[-5]), int(s[-4]), str(s[-3] or ""),
                 str(s[-2]).strip(), str(s[-1] or "")) for s in specs)
  return kind, min(atoms, atoms[::-1])


@contextmanager
def _cylinder_temporary_restraints(mol_id, restraints):
  owned = []
  try:
    existing = {_cylinder_restraint_key(row[0], row[1:3 if row[0] == "bond" else 4])
                for row in (coot.list_extra_restraints_py(mol_id) or []) if row[0] in ("bond", "angle")}
    for kind, specs, target, sigma in restraints:
      key = _cylinder_restraint_key(kind, specs)
      if key in existing:
        continue
      add = coot.add_extra_bond_restraint if kind == "bond" else coot.add_extra_angle_restraint
      index = add(mol_id, *(field for spec in specs for field in spec), target, sigma)
      if index is None or index < 0:
        raise RuntimeError("Could not add a temporary secondary-structure restraint.")
      owned.append([kind, *specs])
      existing.add(key)
    yield len(owned)
  finally:
    errors = []
    for row in reversed(owned):
      try:
        coot.delete_extra_restraint_py(mol_id, row)
      except Exception as error:
        errors.append(str(error))
    if errors:
      raise RuntimeError("Could not remove temporary cylinder restraints: " + "; ".join(errors))


def _start_cylinder_refinement(mol_id, moving, primary, secondary_structure=None):
  if coot.refinement_already_ongoing_p():
    raise ValueError("Accept or reject the current refinement first.")
  if not valid_model_molecule_qm(mol_id) or not valid_map_molecule_qm(imol_refinement_map()):
    raise ValueError("Choose an available model and refinement map.")
  restraints = []
  if secondary_structure:
    required = ("list_extra_restraints_py", "delete_extra_restraint_py", "add_extra_bond_restraint")
    if secondary_structure == "beta strand":
      required += ("add_extra_angle_restraint",)
    if any(not callable(getattr(coot, name, None)) for name in required):
      raise ValueError("This Coot build cannot safely add and remove temporary secondary-structure restraints.")
    restraints = _cylinder_secondary_restraints(_model_coordinate_snapshot(mol_id), primary, secondary_structure)
  previous = (coot.refinement_immediate_replacement_state(), coot.secondary_structure_restraints_type())
  attempted = False
  try:
    coot.set_refinement_immediate_replacement(0)
    # Disable forced global SS. Coot still auto-preserves detected helices in
    # the contact shell; its Python API does not expose a switch for that.
    coot.set_secondary_structure_restraints_type(0)
    with _cylinder_temporary_restraints(mol_id, restraints) if restraints else nullcontext():
      attempted = True
      result = coot.refine_residues_py(mol_id, moving)
      if result is False or result is None or not coot.refinement_already_ongoing_p():
        raise RuntimeError("Coot could not start cylinder refinement; check the selection and dictionaries.")
    # Native refinement owns a copy of the restraints, so the source-model extras
    # can be removed now without losing them from the live intermediate atoms.
    add_status_bar_text(f"Cylinder refinement: {len(primary)} primary + {len(moving)-len(primary)} contact residues"
                        + (f" | {secondary_structure}" if secondary_structure else "")
                        + ". Inspect intermediate atoms and Accept or Reject in Coot.")
    return result
  except Exception:
    if attempted and coot.refinement_already_ongoing_p():
      coot.clear_up_moving_atoms()
    raise
  finally:
    try:
      coot.set_secondary_structure_restraints_type(previous[1])
    finally:
      coot.set_refinement_immediate_replacement(previous[0])


def _pick_cylinder_refinement(secondary_structure=None):
  if not valid_map_molecule_qm(imol_refinement_map()):
    info_dialog("Select a refinement map before cylinder refinement.")
    return
  if coot.refinement_already_ongoing_p():
    info_dialog("Accept or reject the current refinement first.")
    return
  markers = []
  def cleanup():
    while markers:
      try:
        coot.close_generic_object(markers.pop())
      except Exception as error:
        print("Could not clear cylinder endpoint marker:", error)
    coot.graphics_draw()
  def mark(spec, index):
    _mark_register_endpoint(spec, "cyan" if index == 1 else "orange", markers,
                            "Cylinder refinement endpoint")
    add_status_bar_text(_navigation_status_bar_label(_click_spec_imol(spec), _click_spec_chain_id(spec),
      _click_spec_res_no(spec), _click_spec_ins_code(spec)).replace("Current residue:", "Clicked residue:", 1)
      + (" | Click the end residue." if index == 1 else " | Range selected."))
  def picked(first, last):
    try:
      primary = _clicked_residue_specs(first, last)
      if primary:
        return auto_refine(primary, _click_spec_imol(first), secondary_structure)
    except (ValueError, RuntimeError) as error:
      info_dialog(f"Cylinder refinement could not start.\n\n{error}")
    finally:
      cleanup()
  add_status_bar_text("Cylinder refinement: click start (cyan), then end (orange). Tab cancels.")
  return _script_user_defined_click(2, picked, on_pick=mark, cleanup=cleanup)


def cylinder_refine(mol_id=None, chain_id=None, start=None, end=None, secondary_structure=None):
  """Refine a primary range plus local contacts, leaving native intermediate atoms."""
  if all(value is None for value in (mol_id, chain_id, start, end)):
    return _pick_cylinder_refinement(secondary_structure)
  if any(value is None for value in (mol_id, chain_id, start, end)):
    raise ValueError("Supply model, chain, start and end, or use no arguments to pick the range.")
  if mol_id not in model_molecule_list():
    raise ValueError("Select an available coordinate model.")
  try:
    start, end = sorted(int(str(value).strip()) for value in (start, end))
  except (TypeError, ValueError):
    raise ValueError("Range endpoints must be integer residue numbers.") from None
  primary = [tuple(spec) for spec in _tq_chain_residue_specs(mol_id, chain_id) if start <= spec[1] <= end]
  if not primary or not any(s[1] == start for s in primary) or not any(s[1] == end for s in primary):
    raise ValueError("Both range endpoints must exist in the selected chain.")
  return auto_refine(primary, mol_id, secondary_structure)


def cylinder_refine_alpha_helix(mol_id=None, chain_id=None, start=None, end=None):
  return cylinder_refine(mol_id, chain_id, start, end, secondary_structure="alpha helix")


def cylinder_refine_beta_strand(mol_id=None, chain_id=None, start=None, end=None):
  return cylinder_refine(mol_id, chain_id, start, end, secondary_structure="beta strand")
   

#**** "Custom menu item functions ****
#Fits all polymer chains to map
def rigid_fit_all_chains():
  mol_id = _active_molecule_or_status()
  if mol_id is None:
    return None
  with _grouped_model_edit(mol_id):
    for ch_id in chain_ids(mol_id):
      if is_polymer(mol_id, ch_id) == 1:
        rigid_body_refine_by_atom_selection(mol_id, "//%s//" % ch_id)
        accept_regularizement()
  
#Fits active chain to map
def rigid_fit_active_chain():
  chain_context = _active_chain_context_or_status()
  if not chain_context:
    return None
  rigid_body_refine_by_atom_selection(chain_context["mol_id"], "//%s//" % chain_context["chain_id"])
  
#Copies active chain
def copy_active_chain():
  chain_context = _active_chain_context_or_status()
  if not chain_context:
    return None
  new_molecule_by_atom_selection(chain_context["mol_id"], "//%s//" % chain_context["chain_id"])
  
#Cuts active chain
def cut_active_chain():
  chain_context = _active_chain_context_or_status()
  if not chain_context:
    return None
  mol_id = chain_context["mol_id"]
  ch_id = chain_context["chain_id"]
  return _cut_residue_specs_verified(mol_id, _tq_chain_residue_specs(mol_id, ch_id))


def _atom_snapshot_for_specs(mol_id, specs):
  if coot.n_models(mol_id) != 1:
    raise ValueError("This operation requires a single-model molecule.")
  result = {}
  for ch, number, ins in specs:
    atoms = residue_info_py(mol_id, ch, number, ins) or []
    if not atoms:
      raise ValueError(f"Residue {ch}:{number}{ins} is missing.")
    name = residue_name(mol_id, ch, number, ins)
    for atom in atoms:
      key = (ch, number, ins, name, atom[0][0], atom[0][1], repr(atom[1]))
      if key in result:
        raise ValueError("Ambiguous atom identity in selection.")
      result[key] = tuple(atom[2])
  if not result:
    raise ValueError("The selection contains no atoms.")
  return result


def _selection_for_residue_specs(specs):
  return "||".join(_residue_spec_to_cid_compat(spec) for spec in specs)


def _copy_residue_specs_verified(mol_id, specs):
  """Verify an independent, complete copy before any destructive edit."""
  original = _atom_snapshot_for_specs(mol_id, specs)
  existing = set(model_molecule_list())
  copy_function = getattr(coot, "new_molecule_by_residue_specs_py", None)
  copied = (copy_function(mol_id, [list(spec) for spec in specs]) if callable(copy_function)
            else new_molecule_by_atom_selection(mol_id, _selection_for_residue_specs(specs)))
  if copied in existing or not valid_model_molecule_qm(copied):
    raise RuntimeError("Coot could not create an independent copy; original unchanged.")
  try:
    copied_specs = _all_residue_specs_for_colouring(copied)
    if _atom_snapshot_for_specs(copied, copied_specs) != original:
      raise RuntimeError("Copied atoms do not match the selection; original unchanged.")
  except Exception:
    close_molecule(copied)
    raise
  return copied


def _delete_residue_specs(mol_id, specs):
  delete_many = getattr(coot, "delete_residues_py", None)
  if callable(delete_many):
    delete_many(mol_id, [list(spec) for spec in specs])
  else:
    for spec in specs:
      delete_residue(mol_id, *spec)


def _cut_residue_specs_verified(mol_id, specs):
  try:
    copied = _copy_residue_specs_verified(mol_id, specs)
    with _grouped_model_edit(mol_id):
      _delete_residue_specs(mol_id, specs)
    return copied
  except Exception as error:
    info_dialog(f"Cut stopped: {error}")
    return -1

def _rigid_fit_atom_snapshot(mol_id, chain=None, start=None, end=None):
  """Index coordinates by full atom identity, not transient serial indices."""
  if coot.n_models(mol_id) != 1:
    raise ValueError("Fast rigid fitting requires a single-model molecule.")
  snapshot = {}
  for ch, number, ins in _all_residue_specs_for_colouring(mol_id):
    if chain is not None and (ch != chain or not start <= number <= end):
      continue
    name = residue_name(mol_id, ch, number, ins)
    for atom in coot.residue_info_py(mol_id, ch, number, ins) or []:
      key = (ch, number, ins, name, atom[0][0], atom[0][1], repr(atom[1]))
      xyz = tuple(float(v) for v in atom[2])
      if key in snapshot or len(xyz) != 3 or not all(math.isfinite(v) for v in xyz):
        raise ValueError("Fragment has ambiguous atom identities or invalid coordinates.")
      snapshot[key] = xyz
  if not snapshot:
    raise ValueError("No atoms found in the requested range.")
  return snapshot


def fast_rigid_fit(res_start, res_end, ch_id, mol_id):
  """Fit a disposable copy, then update matching coordinates without deleting residues."""
  fragment = None
  previous_replacement = None
  updating_source = False
  try:
    if mol_id not in model_molecule_number_list() or ch_id not in chain_ids(mol_id):
      raise ValueError("The specified chain or molecule does not exist.")
    if coot.imol_refinement_map() not in map_molecule_list():
      raise ValueError("Select a refinement map before rigid fitting.")
    if coot.backup_state(mol_id) != 1:
      raise ValueError("Enable Coot backups for this model before fast rigid fitting.")
    pending = getattr(coot, "get_continue_updating_refinement_atoms_state", None)
    if callable(pending) and pending():
      raise ValueError("Accept or cancel the current refinement before rigid fitting.")
    start, end = sorted((int(res_start), int(res_end)))
    original = _rigid_fit_atom_snapshot(mol_id, ch_id, start, end)
    selection = _atom_selection_for_residue_range(ch_id, start, end)
    existing = set(model_molecule_number_list())
    copied = coot.new_molecule_by_atom_selection(mol_id, selection)
    if copied in existing or not valid_model_molecule_qm(copied):
      raise RuntimeError("Coot could not create an independent fitting fragment.")
    fragment = copied
    coot.set_mol_displayed(fragment, 0)
    if _rigid_fit_atom_snapshot(fragment) != original:
      raise RuntimeError("The copied fragment does not match the requested atoms.")
    previous_replacement = coot.refinement_immediate_replacement_state()
    coot.set_refinement_immediate_replacement(1)
    coot.rigid_body_refine_by_atom_selection(fragment, "//")
    coot.accept_regularizement()
    fitted = _rigid_fit_atom_snapshot(fragment)
    if fitted.keys() != original.keys():
      raise RuntimeError("Rigid fitting changed fragment atom identities.")
    if _rigid_fit_atom_snapshot(mol_id, ch_id, start, end) != original:
      raise RuntimeError("The original region changed during fitting; no replacement applied.")
    with _grouped_model_edit(mol_id):
      updating_source = True
      if not coot.replace_fragment(mol_id, fragment, selection):
        raise RuntimeError("Coot rejected the fragment update.")
      updated = _rigid_fit_atom_snapshot(mol_id, ch_id, start, end)
      if updated.keys() != fitted.keys() or any(
          sum((a-b)**2 for a, b in zip(updated[key], xyz)) > 1e-8 for key, xyz in fitted.items()):
        raise RuntimeError("The updated coordinates do not match the fitted fragment.")
    add_status_bar_text(f"Rigid fitted {ch_id}:{start}-{end}; one Undo restores the original coordinates.")
    return True
  except Exception as error:
    recovery = "\nThe original may have been updated; use Coot Undo to recover." if updating_source else "\nOriginal coordinates were not replaced."
    info_dialog(f"Fast rigid fit failed: {error}{recovery}")
    return False
  finally:
    try:
      if previous_replacement is not None:
        coot.set_refinement_immediate_replacement(previous_replacement)
    finally:
      if fragment is not None and valid_model_molecule_qm(fragment):
        coot.close_molecule(fragment)

#Copy active segment
def copy_active_segment():
  segment_context = _active_segment_context_or_status()
  if not segment_context:
    return None
  return _copy_residue_specs_verified(segment_context["mol_id"], segment_context["segment_specs"])

#Cut active segment
def cut_active_segment():
  segment_context = _active_segment_context_or_status()
  if not segment_context:
    return None
  mol_id = segment_context["mol_id"]
  ch_id = segment_context["chain_id"]
  res_start = segment_context["segment_start"]
  res_end = segment_context["segment_end"]
  return _cut_residue_specs_verified(mol_id, segment_context["segment_specs"])

#Delete active segment
def delete_active_segment():
  segment_context = _active_segment_context_or_status()
  if not segment_context:
    return None
  with _grouped_model_edit(segment_context["mol_id"]):
    _delete_residue_specs(segment_context["mol_id"], segment_context["segment_specs"])


# Random-start rigid fitting and candidate inspection.
_GLOBAL_FIT_STATE = None


def _global_fit_read_grid(path, np):
  """Read Coot's own float CCP4 export, including axis order and non-orthogonal cells."""
  with open(path, "rb") as handle:
    header = handle.read(1024)
  if len(header) != 1024:
    raise ValueError("Incomplete exported map header.")
  for endian in ("<", ">"):
    fields = struct.unpack(endian + "256i", header)
    shape = fields[:3]
    axes = fields[16:19]
    if fields[3] == 2 and all(2 <= n <= 100000 for n in shape) and sorted(axes) == [1, 2, 3]:
      break
  else:
    raise ValueError("Unsupported Coot map export (expected a float CCP4 grid).")
  sampling = np.asarray(fields[7:10], dtype=float)
  cell_values = struct.unpack(endian + "6f", header[40:64])
  a, b, c = cell_values[:3]
  alpha, beta, gamma = map(math.radians, cell_values[3:])
  if min(a, b, c, *sampling) <= 0 or not all(0 < v < math.pi for v in (alpha, beta, gamma)):
    raise ValueError("Invalid exported map cell.")
  cy = c*(math.cos(alpha)-math.cos(beta)*math.cos(gamma))/math.sin(gamma)
  cz_sq = c*c-(c*math.cos(beta))**2-cy*cy
  if cz_sq <= 0:
    raise ValueError("Degenerate exported map cell.")
  cell = np.asarray(((a, b*math.cos(gamma), c*math.cos(beta)),
                     (0., b*math.sin(gamma), cy), (0., 0., math.sqrt(cz_sq))))
  transform = (np.diag(sampling) @ np.linalg.inv(cell))[np.asarray(axes)-1]
  offset = 1024 + fields[23]
  if fields[23] < 0 or os.path.getsize(path) < offset + math.prod(shape)*4:
    raise ValueError("Incomplete exported map data.")
  data = np.memmap(path, dtype=endian + "f4", mode="r", offset=offset, shape=tuple(reversed(shape)))
  periods = sampling[np.asarray(axes)-1].astype(int)
  # Crystallographic exports can span a complete period including duplicate end planes.
  periodic = fields[22] > 1 and all(n >= p for n, p in zip(shape, periods))
  return {"np": np, "data": data, "shape": np.asarray(shape), "transform": transform, "spacegroup": fields[22],
          "starts": np.asarray(fields[4:7], dtype=float), "periodic": periodic, "periods": periods,
          "origin": np.asarray(struct.unpack(endian + "3f", header[196:208]), dtype=float)}


def _global_fit_grid_values(grid, points, gradients=False):
  np = grid["np"]
  q = (points-grid["origin"]) @ grid["transform"].T - grid["starts"]
  valid = np.isfinite(q).all(axis=1)
  if grid["periodic"]:
    q = np.mod(q, grid["periods"])
  else:
    valid &= ((q >= 0) & (q <= grid["shape"]-1)).all(axis=1)
  q = np.where(valid[:, None], q, 0.)
  lower = np.floor(q).astype(int)
  if not grid["periodic"]:
    lower = np.minimum(lower, grid["shape"]-2)
  fx, fy, fz = (q-lower).T
  if grid["periodic"]:
    # Tiny negative coordinates can round to the exact period after modulo.
    lower %= grid["periods"]
  x0, y0, z0 = lower.T
  x1, y1, z1 = (lower+1).T
  if grid["periodic"]:
    x1 = x1 % grid["periods"][0]
    y1 = y1 % grid["periods"][1]
    z1 = z1 % grid["periods"][2]
  data = np.asarray(grid["data"])
  a = data[z0, y0, x0].astype(float, copy=False)
  b = data[z0, y0, x1].astype(float, copy=False)
  c = data[z0, y1, x0].astype(float, copy=False)
  d = data[z0, y1, x1].astype(float, copy=False)
  e = data[z1, y0, x0].astype(float, copy=False)
  f = data[z1, y0, x1].astype(float, copy=False)
  g = data[z1, y1, x0].astype(float, copy=False)
  h = data[z1, y1, x1].astype(float, copy=False)
  # Nested interpolation reuses the same eight corners for values and gradients.
  dx00, dx10, dx01, dx11 = b-a, d-c, f-e, h-g
  x00, x10 = a+fx*dx00, c+fx*dx10
  x01, x11 = e+fx*dx01, g+fx*dx11
  dy0, dy1 = x10-x00, x11-x01
  y_low, y_high = x00+fy*dy0, x01+fy*dy1
  values = y_low+fz*(y_high-y_low)
  values[~valid] = 0.
  if not np.isfinite(values).all():
    raise ValueError("The fitting map contains non-finite density values.")
  if gradients:
    gx0, gx1 = dx00+fy*(dx10-dx00), dx01+fy*(dx11-dx01)
    derivative = np.column_stack((gx0+fz*(gx1-gx0), dy0+fz*(dy1-dy0), y_high-y_low))
    derivative[~valid] = 0.
    derivative = derivative @ grid["transform"]
  else:
    derivative = None
  return values, derivative, valid


def _global_fit_load_grid(map_id, original, radius):
  np = _tq_numpy_module()
  if np is None:
    raise RuntimeError("Global fitting needs NumPy in Coot's Python environment.")
  owner = tempfile.TemporaryDirectory(prefix="coot-global-fit-")
  grid = None
  try:
    path = os.path.join(owner.name, "target.map")
    if not coot.export_map(map_id, path):
      raise RuntimeError("Coot could not export the active map for fitting.")
    grid = _global_fit_read_grid(path, np)
    parameters = getattr(coot, "map_parameters_py", lambda _m: None)(map_id)
    # Coot returns a five-item list even for EM maps, with empty MTZ fields.
    crystal = grid["spacegroup"] > 1 or bool(isinstance(parameters, (list, tuple)) and
                                            len(parameters) >= 3 and all(parameters[:3]))
    if crystal and all(n >= p for n, p in zip(grid["shape"], grid["periods"])):
      grid["periodic"] = True
    points = np.asarray(list(original.values()), dtype=float)
    if not np.isfinite(points).all():
      raise ValueError("The selection has non-finite coordinates.")
    _values, _gradient, inside = _global_fit_grid_values(grid, points)
    if crystal and not grid["periodic"]:
      # CCP4's crystallographic exporter can write only an asymmetric unit.
      # Ask Coot to expand a Cartesian box, preserving the exact map samples.
      centre = points.mean(axis=0)
      extent = float(np.linalg.norm(points-centre, axis=1).max()) + radius + 25.
      sizes = 2*np.ceil(extent*np.linalg.norm(grid["transform"], axis=1)) + 3
      if math.prod(map(int, sizes)) > 128000000:
        raise ValueError("The fitting box would exceed 512 MB. Coarsen the map or choose a smaller selection first.")
      grid["data"]._mmap.close()
      grid = None
      if not coot.export_map_fragment(map_id, *map(float, centre), extent, path):
        raise ValueError("Coot could not expand the fitting map region.")
      grid = _global_fit_read_grid(path, np)
    grid["owner"] = owner
    # Verify header/coordinate interpretation at grid nodes against Coot itself.
    if inside.any():
      available = points[inside]
      probe = available[np.linspace(0, len(available)-1, min(12, len(available)), dtype=int)]
      index = np.rint((probe-grid["origin"]) @ grid["transform"].T-grid["starts"])
    else:
      index = np.asarray(list(product((.25, .75), repeat=3)))*(grid["shape"]-1)
      index = np.rint(index)
    probe = (index+grid["starts"]) @ np.linalg.inv(grid["transform"]).T + grid["origin"]
    values, _gradient, inside = _global_fit_grid_values(grid, probe)
    for p, value, valid in zip(probe, values, inside):
      if not valid:
        raise ValueError("The selection lies outside the exported map bounds; use a boxed map covering the selection.")
      expected = float(coot.density_at_point(map_id, *map(float, p)))
      if not math.isfinite(expected) or abs(expected-value) > .002*max(1., abs(expected)):
        raise ValueError("Exported map coordinates do not match Coot; no search was started.")
    return grid
  except Exception:
    if grid is not None:
      _global_fit_close_grid(grid)
    owner.cleanup()
    raise


def _global_fit_close_grid(grid):
  if grid is not None:
    mapping = getattr(grid.get("data"), "_mmap", None)
    if mapping is not None:
      mapping.close()
    pool = grid.get("worker_pool")
    if pool is not None:
      _global_fit_workers_stop(pool)
    owner = grid.get("owner")
    if owner is not None:
      if pool is not None and not pool["reaped"].is_set():
        import threading
        def cleanup():
          pool["reaped"].wait()
          owner.cleanup()
        threading.Thread(target=cleanup, daemon=True).start()
      else:
        owner.cleanup()


def _global_fit_optimize_iter(snapshot, grid, deadline=None, return_transform=False):
  """Independent rigid gradient ascent on mean trilinearly interpolated atom density."""
  np = grid["np"]
  all_points = np.asarray(snapshot if isinstance(snapshot, np.ndarray) else list(snapshot.values()), dtype=float)
  centre0 = all_points.mean(axis=0)
  centre = centre0.copy()
  # Keep the optimisation cost bounded for large structures; score the complete pose.
  indices = np.linspace(0, len(all_points)-1, min(3000, len(all_points)), dtype=int)
  points = all_points[indices].copy()
  rotation = np.eye(3)
  step = .5
  values, gradients, valid = _global_fit_grid_values(grid, points, gradients=True)
  score = float(values.mean())
  if deadline is None:
    deadline = time.monotonic() + .020
  for _iteration in range(2000):
    if time.monotonic() >= deadline:
      yield _iteration
      deadline = time.monotonic() + .020
    if step < .01:
      break
    proposals = []
    direction = gradients.mean(axis=0)
    norm = np.linalg.norm(grid["transform"] @ direction)
    if norm > 1e-15:
      shift = direction*(step/norm)
      proposed = points+shift
      v, proposed_gradient, _valid = _global_fit_grid_values(grid, proposed, gradients=True)
      proposals.append((float(v.mean()), proposed, shift, None, proposed_gradient))
    offsets = points-centre
    torque = np.cross(offsets, gradients).sum(axis=0)
    torque_norm = np.linalg.norm(torque)
    if torque_norm > 1e-15:
      axis = torque/torque_norm
      lever = float(np.linalg.norm(np.cross(axis, offsets) @ grid["transform"].T, axis=1).max())
      if lever > 1e-12:
        angle = min(step/lever, math.radians(10.))
        x, y, z = axis
        skew = np.asarray(((0., -z, y), (z, 0., -x), (-y, x, 0.)))
        delta = np.eye(3) + math.sin(angle)*skew + (1-math.cos(angle))*(skew @ skew)
        proposed = offsets @ delta.T+centre
        v, proposed_gradient, _valid = _global_fit_grid_values(grid, proposed, gradients=True)
        proposals.append((float(v.mean()), proposed, None, delta, proposed_gradient))
    best = max(proposals, key=lambda p: p[0], default=None)
    if best is None or best[0] <= score + 1e-10*max(1., abs(score)):
      step *= .5
      continue
    score, points, shift, delta, gradients = best
    if shift is not None:
      centre += shift
    else:
      rotation = delta @ rotation
    step = min(.5, step*1.1)
  fitted = (all_points-centre0) @ rotation.T+centre
  values, _gradient, valid = _global_fit_grid_values(grid, fitted)
  if not valid.any():
    return None, None
  if return_transform:
    return (rotation, centre-centre0 @ rotation.T), float(values.mean())
  return dict(zip(snapshot, map(tuple, fitted))), float(values.mean())


def _global_fit_optimize(snapshot, grid):
  """Synchronous numerical entry point; the GUI advances the iterator in short slices."""
  iterator = _global_fit_optimize_iter(snapshot, grid)
  while True:
    try:
      next(iterator)
    except StopIteration as result:
      return result.value


def _global_fit_random_pose(radius, rng):
  """Uniform orientation and uniform-volume displacement of the selection centre."""
  q = [rng.gauss(0.0, 1.0) for _ in range(4)]
  length = math.sqrt(sum(v*v for v in q))
  if length < 1e-12:
    q, length = [1.0, 0.0, 0.0, 0.0], 1.0
  w, x, y, z = [v/length for v in q]
  rotation = ((1-2*(y*y+z*z), 2*(x*y-z*w), 2*(x*z+y*w)),
              (2*(x*y+z*w), 1-2*(x*x+z*z), 2*(y*z-x*w)),
              (2*(x*z-y*w), 2*(y*z+x*w), 1-2*(x*x+y*y)))
  direction = [rng.gauss(0.0, 1.0) for _ in range(3)]
  norm = math.sqrt(sum(v*v for v in direction))
  distance = radius * rng.random()**(1.0/3.0)
  shift = [v*distance/norm if norm > 1e-12 else 0.0 for v in direction]
  return rotation, shift


def _global_fit_random_start(original, radius, rng):
  """Uniform orientations and uniform-volume translations about the original centroid."""
  rotation, shift = _global_fit_random_pose(radius, rng)
  centre = [sum(p[i] for p in original.values())/len(original) for i in range(3)]
  return {key: tuple(centre[i] + shift[i] + sum(rotation[i][j]*(p[j]-centre[j]) for j in range(3))
                     for i in range(3)) for key, p in original.items()}


def _global_fit_worker_count():
  count = getattr(os, "process_cpu_count", os.cpu_count)() or 1
  if callable(getattr(os, "sched_getaffinity", None)):
    try:
      count = min(count, len(os.sched_getaffinity(0)))
    except OSError:
      pass
  return max(1, min(8, count-1))


def _global_fit_worker_python():
  version = f"python{sys.version_info.major}.{sys.version_info.minor}"
  candidates = [os.path.join(prefix, "bin", name) for prefix in (sys.prefix, sys.base_prefix)
                for name in (version, "python3")]
  candidates.extend((os.path.join(sys.prefix, "python.exe"), shutil.which(version), shutil.which("python3")))
  for path in candidates:
    if path and os.path.isfile(path) and os.access(path, os.X_OK):
      return path
  raise RuntimeError("No standalone Python interpreter found; using one fitting worker.")


def _global_fit_worker_main():
  """Pure numerical subprocess entry point. No Coot or GTK imports or calls."""
  import json
  import os
  import sys
  import traceback
  import numpy as np
  try:
    with open(sys.argv[1]) as handle:
      config = json.load(handle)
    data = np.memmap(config.pop("path"), dtype=config.pop("dtype"), mode="r",
                     offset=config.pop("offset"), shape=tuple(config.pop("data_shape")))
    original = np.load(config.pop("original"), mmap_mode="r")
    grid = dict(config, np=np, data=data)
    for key in ("shape", "transform", "starts", "periods", "origin"):
      grid[key] = np.asarray(grid[key])
    centre = original.mean(axis=0)
    parent = os.getppid()
    print(json.dumps({"ready": True}), flush=True)
    for line in sys.stdin:
      job = json.loads(line)
      rotation = np.asarray(job["rotation"])
      shift = np.asarray(job["shift"])
      translation = centre+shift-centre @ rotation.T
      points = original @ rotation.T+translation
      iterator = _global_fit_optimize_iter(points, grid, return_transform=True)
      last_progress = 0.
      while True:
        if os.getppid() != parent:
          return
        try:
          iteration = next(iterator)
          if time.monotonic()-last_progress >= .1:
            print(json.dumps({"trial": job["trial"], "iteration": iteration}), flush=True)
            last_progress = time.monotonic()
        except StopIteration as result:
          transform, score = result.value
          if transform is None:
            response = {"trial": job["trial"], "score": None}
          else:
            local_rotation, local_translation = transform
            response = {"trial": job["trial"], "score": score,
                        "rotation": (local_rotation @ rotation).tolist(),
                        "translation": (translation @ local_rotation.T+local_translation).tolist()}
          print(json.dumps(response, allow_nan=False), flush=True)
          break
  except BrokenPipeError:
    return
  except Exception as error:
    traceback.print_exc(file=sys.stderr)
    print(json.dumps({"fatal": str(error)}), flush=True)


def _global_fit_worker_source():
  import inspect
  functions = (_global_fit_grid_values, _global_fit_optimize_iter, _global_fit_worker_main)
  return "import math, time\n" + "\n".join(inspect.getsource(f) for f in functions) + "\n_global_fit_worker_main()\n"


def _global_fit_workers_stop(pool):
  if pool.get("stopped"):
    return
  import threading
  pool["stopped"] = True
  for worker in pool["workers"]:
    process = worker["process"]
    try:
      process.stdin.close()
    except (OSError, ValueError):
      pass
    try:
      if process.poll() is None:
        process.terminate()
    except OSError:
      pass
  def reap():
    try:
      for worker in pool["workers"]:
        process = worker["process"]
        try:
          process.wait(timeout=1.)
        except Exception:
          try:
            process.kill()
            process.wait()
          except OSError:
            pass
        try:
          process.stdout.close()
        except (OSError, ValueError):
          pass
    finally:
      pool["reaped"].set()
  threading.Thread(target=reap, daemon=True).start()


def _global_fit_workers_read(pool):
  for worker in pool["workers"]:
    process = worker["process"]
    try:
      data = os.read(process.stdout.fileno(), 65536)
    except BlockingIOError:
      data = None
    if data:
      worker["buffer"] += data
      while b"\n" in worker["buffer"]:
        line, worker["buffer"] = worker["buffer"].split(b"\n", 1)
        message = json.loads(line)
        if "fatal" in message:
          raise RuntimeError("Fitting worker failed: " + message["fatal"])
        pool["events"].append((worker, message))
    elif data == b"" or process.poll() is not None:
      raise RuntimeError("A fitting worker exited unexpectedly. Reset and restart the search.")


def _global_fit_workers_prepare(state):
  import subprocess
  import threading
  from collections import deque
  grid = state["grid"]
  np = grid["np"]
  directory = grid["owner"].name
  original = np.asarray(list(state["original"].values()), dtype=float)
  original_path = os.path.join(directory, "original.npy")
  np.save(original_path, original)
  config = {key: grid[key].tolist() for key in ("shape", "transform", "starts", "periods", "origin")}
  config.update(periodic=grid["periodic"], path=str(grid["data"].filename),
                dtype=grid["data"].dtype.str, offset=int(grid["data"].offset),
                data_shape=list(grid["data"].shape), original=original_path)
  config_path = os.path.join(directory, "workers.json")
  with open(config_path, "w") as handle:
    json.dump(config, handle)
  source, python = _global_fit_worker_source(), _global_fit_worker_python()
  env = os.environ.copy()
  for key in ("PYTHONHOME", "PYTHONPATH", "PYTHONSTARTUP", "DYLD_INSERT_LIBRARIES"):
    env.pop(key, None)
  for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS"):
    env[key] = "1"
  pool = dict(workers=[], events=deque(), reaped=threading.Event(), next=0,
              original=original, centre=original.mean(axis=0), keys=tuple(state["original"]))
  offsets = original-pool["centre"]
  pool["covariance"] = offsets.T @ offsets/len(offsets)
  grid["worker_pool"] = state["pool"] = pool
  for index in range(state["workers"]):
    state["message"] = f"Starting fitting workers: {index+1}/{state['workers']}..."
    yield
    with open(os.path.join(directory, f"worker-{index}.log"), "wb") as log:
      process = subprocess.Popen([python, "-E", "-u", "-c", source, config_path],
                                 stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=log, env=env, cwd=directory)
    worker = dict(process=process, buffer=b"", ready=False, job=None)
    pool["workers"].append(worker)
    os.set_blocking(process.stdout.fileno(), False)
  deadline = time.monotonic()+30.
  while not all(w["ready"] for w in pool["workers"]):
    _global_fit_workers_read(pool)
    while pool["events"]:
      worker, message = pool["events"].popleft()
      if message.get("ready") is not True:
        raise RuntimeError("Unexpected fitting-worker startup response.")
      worker["ready"] = True
    if time.monotonic() > deadline:
      raise RuntimeError("Fitting workers did not start in time; using one worker.")
    ready = sum(w["ready"] for w in pool["workers"])
    state["message"] = f"Starting fitting workers: {ready}/{state['workers']} ready..."
    yield


def _global_fit_workers_step(state):
  pool, np = state["pool"], state["grid"]["np"]
  _global_fit_workers_read(pool)
  if pool["events"]:
    worker, result = pool["events"].popleft()
    if result.get("trial") != worker["job"]:
      raise RuntimeError("Unexpected fitting-worker trial identifier.")
    if "iteration" in result:
      worker["iteration"] = state["iteration"] = int(result["iteration"])
    else:
      score = result["score"]
      if score is None:
        state["failed"] += 1
      else:
        rotation = np.asarray(result["rotation"], dtype=float)
        translation = np.asarray(result["translation"], dtype=float)
        if (rotation.shape != (3, 3) or translation.shape != (3,) or
            not np.isfinite(rotation).all() or not np.isfinite(translation).all() or not math.isfinite(score) or
            not np.allclose(rotation @ rotation.T, np.eye(3), atol=1e-6, rtol=0.) or
            abs(float(np.linalg.det(rotation))-1.) > 1e-6):
          raise RuntimeError("Invalid rigid transform from a fitting worker.")
        _global_fit_insert_transform(state, rotation, translation, float(score))
      worker["job"] = None
      worker["iteration"] = 0
      state["done"] += 1
      state["iteration"] = 0
  for worker in pool["workers"]:
    if worker["job"] is None and pool["next"] < state["trials"]:
      trial = pool["next"]
      if trial == 0:
        rotation, shift = ((1., 0., 0.), (0., 1., 0.), (0., 0., 1.)), (0., 0., 0.)
      else:
        rotation, shift = _global_fit_random_pose(state["radius"], state["rng"])
      process = worker["process"]
      process.stdin.write((json.dumps(dict(trial=trial, rotation=rotation, shift=shift))+"\n").encode())
      process.stdin.flush()
      worker["job"] = trial
      worker["iteration"] = 0
      pool["next"] += 1
  state["parallel_fraction"] = sum(min(.99, w.get("iteration", 0)/2000)
                                   for w in pool["workers"] if w["job"] is not None)
  if state["done"] >= state["trials"]:
    _global_fit_workers_stop(pool)
    state.update(phase="score", scored=0, iteration=0,
                 message=f"Finished {state['done']} trials; {state['failed']} without a fit.")


def _global_fit_rmsd(first, second, cutoff=None):
  # No superposition: two copies with the same shape at different placements are distinct.
  total = 0.0
  limit = float("inf") if cutoff is None else len(first)*cutoff*cutoff
  for key, p in first.items():
    total += sum((p[i]-second[key][i])**2 for i in range(3))
    if total > limit:
      return cutoff+1.0
  return math.sqrt(total/len(first))


def _global_fit_insert_transform(state, rotation, translation, score):
  """Cluster rigid placements exactly from second moments, without rebuilding all atoms."""
  pool, np = state["pool"], state["grid"]["np"]
  centre = pool["centre"] @ rotation.T+translation
  match = None
  for candidate in state["candidates"]:
    delta = rotation-candidate["rotation"]
    mse = float(np.sum((delta @ pool["covariance"])*delta) + np.sum((centre-candidate["centre"])**2))
    if max(0., mse) < 1.:
      match = candidate
      break
  pose = dict(rotation=rotation, translation=translation, centre=centre,
              coordinates=None, score=score, cc=None)
  if match is None:
    state["candidates"].append(dict(pose, hits=1))
  else:
    match["hits"] += 1
    if score > match["score"]:
      match.update(pose)
  state["candidates"].sort(key=lambda c: c["score"], reverse=True)
  del state["candidates"][20:]


def _global_fit_candidate_coordinates(state, candidate):
  if candidate["coordinates"] is None:
    pool = state["pool"]
    points = pool["original"] @ candidate["rotation"].T+candidate["translation"]
    candidate["coordinates"] = dict(zip(pool["keys"], map(tuple, points)))
  return candidate["coordinates"]


def _global_fit_write_coordinates(model, snapshot):
  attributes = []
  for key, xyz in snapshot.items():
    if len(xyz) != 3 or not all(math.isfinite(v) for v in xyz):
      raise ValueError("Invalid fitting coordinates.")
    chain, number, ins, _name, atom, alt, _properties = key
    for axis, value in zip(("x", "y", "z"), xyz):
      attributes.append([model, chain, number, ins, atom, alt, axis, float(value)])
  # This API returns zero even on success; verify coordinates at acceptance instead.
  coot.set_atom_attributes_py(attributes)


def _global_fit_insert_candidate(state, snapshot, score):
  if not math.isfinite(score):
    raise RuntimeError("Map density is unavailable for this placement.")
  match = next((c for c in state["candidates"] if _global_fit_rmsd(snapshot, c["coordinates"], cutoff=1.0) < 1.0), None)
  if match is not None:
    match["hits"] += 1
    if score > match["score"]:
      match.update(coordinates=snapshot, score=score, cc=None)
  else:
    state["candidates"].append({"coordinates": snapshot, "score": score, "cc": None, "hits": 1})
  state["candidates"].sort(key=lambda c: c["score"], reverse=True)
  del state["candidates"][20:]


def _global_fit_preview_owned(state, key="preview"):
  model = state.get(key)
  return (model is not None and valid_model_molecule_qm(model) and
          coot.molecule_name(model) == state.get(key + "_name"))


def _global_fit_release(state):
  for key in ("preparation", "trial"):
    iterator = state.pop(key, None)
    if iterator is not None:
      iterator.close()
  for key in ("preview", "scorer"):
    if _global_fit_preview_owned(state, key):
      coot.close_molecule(state[key])
    state[key] = None
  _global_fit_close_grid(state.get("grid"))
  state["grid"] = None
  state["pool"] = None


def _global_fit_prepare(state):
  fitting_map = state["map"]
  blurred = None
  blur_name = None
  try:
    if state["blur_b"]:
      state["message"] = f"Preparing fitting map: B-factor blur {state['blur_b']:g} A^2..."
      yield
      existing = set(map_molecule_list())
      model = coot.sharpen_blur_map(fitting_map, state["blur_b"])
      if model in existing or not coot.is_valid_map_molecule(model):
        raise RuntimeError("Could not create a temporary blurred map; target map unchanged.")
      blurred = model
      blur_name = f"Global fit blurred map (temporary) #{model}"
      try:
        coot.set_molecule_name(model, blur_name)
      except Exception:
        coot.close_molecule(model)
        blurred = None
        raise
      coot.set_map_displayed(model, 0)
      fitting_map = model
    state["message"] = "Caching the fitting map..."
    yield
    state["grid"] = _global_fit_load_grid(fitting_map, state["original"], state["radius"])
  finally:
    if (blurred is not None and coot.is_valid_map_molecule(blurred) and
        coot.molecule_name(blurred) == blur_name):
      coot.close_molecule(blurred)
  state["message"] = "Preparing the preview..."
  yield
  for key in ("preview", "scorer"):
    model = _copy_residue_specs_verified(state["source"], state["specs"])
    name = f"Global fit {key} (temporary) #{model}"
    try:
      coot.set_molecule_name(model, name)
    except Exception:
      coot.close_molecule(model)
      raise
    state[key], state[key + "_name"] = model, name
    coot.set_mol_displayed(model, 0)
    coot.set_mol_active(model, 0)
    coot.turn_off_backup(model)
  if not state.get("hidden_original") and coot.mol_is_displayed(state["source"]):
    coot.set_mol_displayed(state["source"], 0)
    state["hidden_original"] = True
  _global_fit_choose(state)
  if state["workers"] > 1:
    try:
      yield from _global_fit_workers_prepare(state)
    except Exception as error:
      pool = state.pop("pool", None)
      if pool is not None:
        _global_fit_workers_stop(pool)
      state["workers"] = 1
      state["worker_note"] = f"Parallel workers unavailable: {error} Using one worker."
      state["message"] = state["worker_note"]
      print("Global fit:", state["message"])
      yield


def _global_fit_schedule(state):
  if state["running"] and not state["closed"]:
    # One-shot, low-priority work gives GTK time to paint between fitting slices.
    glib = _coot_gui_repository_module("GLib")
    state["timer"] = glib.timeout_add(16, _global_fit_dispatch, state,
                                    priority=glib.PRIORITY_DEFAULT_IDLE)


def _global_fit_dispatch(state):
  state["timer"] = None
  if state["closed"] or not state["running"]:
    return False
  state["work_deadline"] = time.monotonic() + .020
  try:
    while _global_fit_step(state, update=False):
      # Batch only completed short trials, within the same 20 ms work budget.
      if (state["phase"] != "search" or state.get("trial") is not None or
          (state.get("pool") is not None and not state["pool"]["events"]) or
          time.monotonic() >= state["work_deadline"]):
        break
  finally:
    state.pop("work_deadline", None)
    state["update"]()
  if state["running"]:
    _global_fit_schedule(state)
  return False


def _global_fit_complete(state):
  if _global_fit_preview_owned(state, "scorer"):
    coot.close_molecule(state["scorer"])
  state["scorer"] = None
  _global_fit_close_grid(state.get("grid"))
  state["grid"] = None
  state["pool"] = None
  state.update(running=False, phase="ready")
  if state["candidates"]:
    _global_fit_choose(state, state["candidates"][0])


def _global_fit_choose(state, candidate=None):
  if state["closed"] or state.get("preview") is None:
    return
  preview = state["preview"]
  if not _global_fit_preview_owned(state):
    raise RuntimeError("The fitting preview was closed. Close and restart the search.")
  _global_fit_write_coordinates(preview, state["original"] if candidate is None else candidate["coordinates"])
  coot.set_mol_displayed(preview, 1)
  state["selected"] = candidate


def _global_fit_stop(state, message=None):
  if not state["running"] or state.get("phase") == "score":
    return
  if state["timer"] is not None:
    _coot_gui_repository_module("GLib").source_remove(state["timer"])
    state["timer"] = None
  state["message"] = message or f"Stopped after {state['done']}/{state['trials']} trials."
  if state.get("pool") is not None:
    _global_fit_workers_stop(state["pool"])
  iterator = state.pop("trial", None)
  if iterator is not None:
    iterator.close()
  if state.get("phase") == "prepare":
    _global_fit_release(state)
    state.update(running=False, phase="ready")
    if state.get("hidden_original") and valid_model_molecule_qm(state["source"]):
      coot.set_mol_displayed(state["source"], 1)
      state["hidden_original"] = False
  else:
    state.update(phase="score", scored=0)
    _global_fit_schedule(state)
  state["update"]()


def _global_fit_step(state, update=True):
  if state["closed"] or not state["running"]:
    return False
  try:
    if not valid_model_molecule_qm(state["source"]):
      raise RuntimeError("The source model was closed.")
    if state["phase"] != "prepare" and not _global_fit_preview_owned(state):
      raise RuntimeError("The source or preview was closed.")
    if not coot.is_valid_map_molecule(state["map"]):
      raise RuntimeError("The fitting map was closed.")
    if coot.refinement_already_ongoing_p():
      raise RuntimeError("Another refinement started. Accept or cancel it before restarting the search.")
    if state["phase"] == "prepare":
      try:
        next(state["preparation"])
      except StopIteration:
        state.pop("preparation")
        state.update(phase="search", message=state.get("worker_note", ""))
    elif state["phase"] == "score":
      if state["scored"] < len(state["candidates"]):
        if not _global_fit_preview_owned(state, "scorer"):
          raise RuntimeError("The temporary scoring model was closed.")
        candidate = state["candidates"][state["scored"]]
        _global_fit_write_coordinates(state["scorer"], _global_fit_candidate_coordinates(state, candidate))
        # Score once against the original target, never against the fitting blur.
        try:
          cc = float(coot.map_to_model_correlation_py(state["scorer"],
                     [list(s) for s in state["specs"]], [], 0, state["map"]))
          if math.isfinite(cc) and -1.0001 <= cc <= 1.0001:
            candidate["cc"] = cc
        except Exception as error:
          print("Global fit: additional CC unavailable:", error)
        state["scored"] += 1
      else:
        _global_fit_complete(state)
    elif state.get("pool") is not None:
      _global_fit_workers_step(state)
    else:
      if state.get("trial") is None:
        coordinates = (state["original"] if state["done"] == 0 else
                       _global_fit_random_start(state["original"], state["radius"], state["rng"]))
        state["trial"] = _global_fit_optimize_iter(coordinates, state["grid"], deadline=state.get("work_deadline"))
      try:
        state["iteration"] = next(state["trial"])
      except StopIteration as result:
        state.pop("trial")
        snapshot, score = result.value
        if snapshot is not None and snapshot.keys() != state["original"].keys():
          raise RuntimeError("Rigid fitting changed atom identities.")
        if snapshot is not None and math.isfinite(score):
          _global_fit_insert_candidate(state, snapshot, score)
        else:
          state["failed"] += 1
        state["done"] += 1
        state["iteration"] = 0
        if state["done"] >= state["trials"]:
          state.update(phase="score", scored=0,
                       message=f"Finished {state['done']} trials; {state['failed']} without a fit.")
  except Exception as error:
    state.update(message=str(error), running=False, phase="error", selected=None)
    _global_fit_release(state)
    if state.get("hidden_original") and valid_model_molecule_qm(state["source"]):
      coot.set_mol_displayed(state["source"], 1)
      state["hidden_original"] = False
    print("Global fit stopped:", error)
  if update:
    state["update"]()
  return state["running"]


def _global_fit_move_to_center(state):
  if state["closed"] or state["running"] or state.get("phase", "idle") != "idle":
    raise ValueError("Reset the search before moving the molecule.")
  source = state["source"]
  if not valid_model_molecule_qm(source):
    raise ValueError("The source model was closed.")
  if coot.refinement_already_ongoing_p():
    raise ValueError("Accept or cancel the current intermediate atoms first.")
  centre = [float(coot.molecule_centre_internal(source, i)) for i in range(3)]
  target = [float(coot.rotation_centre_position(i)) for i in range(3)]
  if not all(math.isfinite(v) for v in centre + target):
    raise ValueError("Cannot determine the molecule or view centre.")
  shift = [t-c for t, c in zip(target, centre)]
  if any(abs(v) > 1e-6 for v in shift):
    with _grouped_model_edit(source):
      coot.translate_molecule_by(source, *shift)
  state["message"] = f"Molecule #{source} moved to screen center. Search will start from this position."
  state["update"]()


def _global_fit_reset(state):
  if state["closed"] or state["running"]:
    raise ValueError("Stop the search and wait for scoring to finish before resetting.")
  if state["timer"] is not None:
    _coot_gui_repository_module("GLib").source_remove(state["timer"])
    state["timer"] = None
  _global_fit_release(state)
  if state.get("hidden_original") and valid_model_molecule_qm(state["source"]):
    coot.set_mol_displayed(state["source"], 1)
  state.update(phase="idle", candidates=[], selected=None, hidden_original=False, done=0,
               scored=0, iteration=0, message="Search cleared. Reposition the molecule or start another search.")
  state["update"]()


def _global_fit_start(state, scope, trials, radius, blur_b=0):
  import random
  if state["running"] or state["closed"]:
    return
  trials, radius, blur_b = int(trials), float(radius), float(blur_b)
  if not 1 <= trials <= 10000 or not math.isfinite(radius) or not 0 <= radius <= 1000:
    raise ValueError("Use 1-10000 trials and a radius between 0 and 1000 A.")
  if not math.isfinite(blur_b) or not 0 <= blur_b <= 100000:
    raise ValueError("Use a blur B-factor between 0 and 100000 A^2 (0 means no blur).")
  if blur_b and not callable(getattr(coot, "sharpen_blur_map", None)):
    raise ValueError("This Coot build cannot create a temporary blurred map.")
  if coot.refinement_already_ongoing_p():
    raise ValueError("Accept or cancel the current intermediate atoms first.")
  if not coot.is_valid_map_molecule(state["map"]):
    raise ValueError("The chosen map is no longer available.")
  source, active = state["source"], state["active"]
  if scope not in ("molecule", "chain", "residue"):
    raise ValueError("Choose molecule, chain or residue.")
  specs = [tuple(s) for s in _all_residue_specs_for_colouring(source)
           if scope == "molecule" or (s[0] == active[0] and (scope == "chain" or tuple(s) == active))]
  original = _atom_snapshot_for_specs(source, specs)
  _global_fit_release(state)
  state.update(preview=None, grid=None, original=original, specs=specs, scope=scope,
               trials=trials, radius=radius, rng=random.Random(), done=0, failed=0,
               candidates=[], selected=None, message="Preparing the search...", running=True,
               phase="prepare", blur_b=blur_b, scored=0, iteration=0,
               workers=min(trials, _global_fit_worker_count()), pool=None, worker_note="")
  state["preparation"] = _global_fit_prepare(state)
  _global_fit_schedule(state)
  state["update"]()


def _global_fit_finish(state, accept=False):
  global _GLOBAL_FIT_STATE
  if state["closed"]:
    return
  if accept:
    if state["running"] or state["selected"] is None:
      raise ValueError("Stop the search and select a solution first.")
    if coot.refinement_already_ongoing_p():
      raise ValueError("Accept or cancel other intermediate atoms first.")
    source = state["source"]
    if not valid_model_molecule_qm(source) or _atom_snapshot_for_specs(source, state["specs"]) != state["original"]:
      raise ValueError("The selected source atoms changed. Close and restart; no fit was applied.")
    current_specs = {tuple(s) for s in _all_residue_specs_for_colouring(source)
                     if state["scope"] == "molecule" or
                     (s[0] == state["active"][0] and (state["scope"] == "chain" or tuple(s) == state["active"]))}
    if current_specs != set(state["specs"]):
      raise ValueError("The source selection changed. Close and restart; no fit was applied.")
    if coot.backup_state(source) != 1:
      raise ValueError("Enable Coot backups for this model so accepting the fit can be undone.")
    snapshot = state["selected"]["coordinates"]
    with _grouped_model_edit(source):
      _global_fit_write_coordinates(source, snapshot)
    actual = _atom_snapshot_for_specs(source, state["specs"])
    if actual.keys() != snapshot.keys() or _global_fit_rmsd(actual, snapshot) > 0.002:
      raise RuntimeError("The coordinate update could not be verified. Use Coot Undo to restore the original.")
    add_status_bar_text(f"Global fit accepted (mean density {state['selected']['score']:.4g}); one Undo restores the original.")
  if state["timer"] is not None:
    _coot_gui_repository_module("GLib").source_remove(state["timer"])
    state["timer"] = None
  state.update(closed=True, running=False)
  _global_fit_release(state)
  if state.get("hidden_original") and valid_model_molecule_qm(state["source"]):
    coot.set_mol_displayed(state["source"], 1)
  state["window"].destroy()
  if _GLOBAL_FIT_STATE is state:
    _GLOBAL_FIT_STATE = None


def global_fit_map(scope="molecule"):
  """Random-start rigid fitting, with an unchanged source and a non-modal solution browser."""
  global _GLOBAL_FIT_STATE
  if _GLOBAL_FIT_STATE is not None:
    _GLOBAL_FIT_STATE["window"].present()
    return _GLOBAL_FIT_STATE
  active = _active_residue_or_status()
  map_id = _scrollable_map_or_status()
  if not active or map_id is None:
    return None
  required = ("export_map", "density_at_point", "map_to_model_correlation_py", "set_atom_attributes_py")
  if any(not callable(getattr(coot, name, None)) for name in required):
    info_dialog("This Coot build does not provide the fitting/coordinate APIs needed for global fitting.")
    return None
  window = Gtk.Window(title="Global fit to map")
  window.set_default_size(540, 440)
  box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=8)
  for side in ("start", "end", "top", "bottom"):
    getattr(box, "set_margin_" + side)(12)
  model_name = str(coot.molecule_name(active[0]) or "Unnamed model")
  map_name = str(coot.molecule_name(map_id) or "Unnamed map")
  label = Gtk.Label(label=f"Model #{active[0]}: {os.path.basename(model_name)} | Map #{map_id}: {os.path.basename(map_name)}\n"
                          f"Active residue {active[1]}:{active[2]}{active[3]}")
  label.set_xalign(0)
  label.set_wrap(True)
  label.set_wrap_mode(_coot_gui_repository_module("Pango").WrapMode.WORD_CHAR)
  label.set_max_width_chars(72)
  label.set_tooltip_text(f"Model #{active[0]}: {model_name}\nMap #{map_id}: {map_name}")
  box.append(label)
  controls = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
  selection = Gtk.ComboBoxText()
  for value in ("molecule", "chain", "residue"):
    selection.append(value, value.capitalize())
  selection.set_active_id(scope if scope in ("molecule", "chain", "residue") else "molecule")
  controls.append(selection)
  fields = []
  for text, value, width in (("Trials", "300", 5), ("Radius (A)", "5", 5)):
    controls.append(Gtk.Label(label=text))
    entry = Gtk.Entry()
    entry.set_text(value)
    entry.set_width_chars(width)
    entry.set_max_width_chars(width)
    controls.append(entry)
    fields.append(entry)
  controls.set_tooltip_text("Random orientations; starting centres within this radius of the original centre. Local fitting can move farther. The first trial fits the original placement.")
  box.append(controls)
  blur_controls = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
  blur = Gtk.CheckButton(label="Blur fitting map")
  blur_b = Gtk.Entry()
  blur_b.set_text("100")
  blur_b.set_width_chars(6)
  blur_b.set_max_width_chars(6)
  blur_b.set_sensitive(False)
  blur.connect("toggled", lambda *_: blur_b.set_sensitive(blur.get_active()))
  blur_controls.append(blur)
  blur_controls.append(Gtk.Label(label="B (A^2)"))
  blur_controls.append(blur_b)
  blur_controls.set_tooltip_text("Optional positive B-factor blur for fitting only. CC is always calculated against the original target map; temporary maps are removed automatically.")
  box.append(blur_controls)
  positioning = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=6)
  move_button = Gtk.Button(label="Move molecule to screen center")
  move_button.set_tooltip_text(f"Translate the whole molecule #{active[0]} to the current view center, preserving its orientation. Available before searching; use Reset search to try another starting point.")
  reset_button = Gtk.Button(label="Reset search")
  reset_button.set_tooltip_text("Discard the search hits without applying a fit, so you can reposition the molecule and search again.")
  positioning.append(move_button)
  positioning.append(reset_button)
  box.append(positioning)
  progress = Gtk.ProgressBar(show_text=True)
  box.append(progress)
  rows = Gtk.ListBox()
  rows.set_selection_mode(Gtk.SelectionMode.SINGLE)
  rows.set_focusable(True)
  groups = [Gtk.SizeGroup(mode=Gtk.SizeGroupMode.HORIZONTAL) for _ in range(4)]
  def table_row(values):
    row = Gtk.ListBoxRow()
    grid = Gtk.Grid(column_spacing=12)
    cells = []
    for index, value in enumerate(values):
      cell = Gtk.Label(label=value)
      cell.set_xalign(1)
      cell.set_width_chars((4, 13, 11, 6)[index])
      cell.add_css_class("monospace")
      groups[index].add_widget(cell)
      grid.attach(cell, index, 0, 1, 1)
      cells.append(cell)
    row.set_child(grid)
    row.global_fit_cells = cells
    return row
  header = Gtk.ListBox()
  header.set_selection_mode(Gtk.SelectionMode.NONE)
  heading = table_row(("Rank", "Mean density", "CC [-1,1]", "Hits"))
  heading.set_selectable(False)
  heading.set_activatable(False)
  heading.global_fit_cells[1].set_tooltip_text("Mean atom density on the fitting map (blurred if enabled). Solutions are ranked by this value.")
  heading.global_fit_cells[2].set_tooltip_text("Pearson correlation against the original target map, from -1 to +1. A dash means unavailable.")
  heading.global_fit_cells[3].set_tooltip_text("Number of trials assigned to this solution cluster; not a confidence estimate.")
  header.append(heading)
  box.append(header)
  scroll = Gtk.ScrolledWindow()
  scroll.set_policy(Gtk.PolicyType.NEVER, Gtk.PolicyType.AUTOMATIC)
  scroll.set_vexpand(True)
  scroll.set_child(rows)
  box.append(scroll)
  status = Gtk.Label(label="Rigid-body search; no flexible refinement. Nearby model atoms are not used to reject fits.")
  status.set_xalign(0)
  status.set_wrap(True)
  box.append(status)
  buttons = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=6)
  actions = {}
  for key, text in (("start", "Search"), ("stop", "Stop"), ("previous", "Previous"), ("next", "Next"),
                    ("original", "Original"), ("accept", "Accept"), ("close", "Close")):
    button = Gtk.Button(label=text)
    buttons.append(button)
    actions[key] = button
  actions.update(move=move_button, reset=reset_button)
  box.append(buttons)
  state = {"window": window, "source": active[0], "active": tuple(active[1:4]), "map": map_id,
           "closed": False, "running": False, "timer": None, "grid": None, "preview": None,
           "selected": None, "candidates": [], "hidden_original": False, "done": 0, "trials": 300,
           "message": "", "rendering": False, "phase": "idle", "rows_signature": None}
  def attempt(function, *args):
    try:
      function(*args)
    except Exception as error:
      status.set_text(str(error))
      print("Global fit:", error)
  def update():
    state["rendering"] = True
    try:
      displayed_candidates = [] if state["running"] else state["candidates"]
      signature = tuple((id(c), c["score"], c["cc"], c["hits"]) for c in displayed_candidates)
      if signature != state["rows_signature"]:
        while rows.get_first_child() is not None:
          rows.remove(rows.get_first_child())
        for index, candidate in enumerate(displayed_candidates):
          cc_text = f"{candidate['cc']:.3f}" if candidate["cc"] is not None else "-"
          row = table_row((str(index+1), f"{candidate['score']:.4g}", cc_text, str(candidate["hits"])))
          row.global_fit_candidate = candidate
          rows.append(row)
        state["rows_signature"] = signature
      selected_row = next((rows.get_row_at_index(i) for i, c in enumerate(displayed_candidates)
                           if c is state["selected"]), None)
      rows.select_row(selected_row)
      running = state["running"]
      available = bool(state["candidates"]) and state.get("preview") is not None
      for control in (selection, *fields, blur):
        control.set_sensitive(not running)
      blur_b.set_sensitive(not running and blur.get_active())
      for key in ("previous", "next"):
        actions[key].set_sensitive(available and not running)
      actions["accept"].set_sensitive(state["selected"] is not None and not running)
      actions["original"].set_sensitive(state["preview"] is not None and not running)
      actions["start"].set_sensitive(not running)
      actions["stop"].set_sensitive(running and state["phase"] != "score")
      rows.set_sensitive(available and not running)
      phase = state["phase"]
      actions["move"].set_sensitive(not running and phase == "idle" and valid_model_molecule_qm(state["source"]))
      actions["reset"].set_sensitive(not running and phase != "idle")
      partial = (state.get("parallel_fraction", 0.) if state.get("pool") is not None else
                 min(.99, state.get("iteration", 0)/2000))
      fraction = .9*(state["done"]+partial)/state["trials"]
      text = f"{state['done']}/{state['trials']} trials | {len(state['candidates'])} distinct solutions"
      if phase == "search":
        text += f" | {state['workers']} worker{'s' if state['workers'] != 1 else ''}"
      if phase == "prepare":
        fraction, text = 0., "Preparing the search..."
      elif phase == "search" and state.get("iteration"):
        text += f" | step {state['iteration']}"
      elif phase == "score":
        fraction = .9 + .1*state["scored"]/max(1, len(state["candidates"]))
        text = f"Scoring solutions: {state['scored']}/{len(state['candidates'])} (original map)"
      elif phase == "ready":
        fraction = 1. if state["done"] else 0.
      progress.set_fraction(fraction)
      progress.set_text(text)
      status.set_text(state["message"] or "Ranked by mean fitting-map density. CC: original map, -1 to +1.")
    finally:
      state["rendering"] = False
  def choose(_rows, row):
    if row is not None and not state["rendering"]:
      attempt(_global_fit_choose, state, row.global_fit_candidate)
      update()
  def navigate(direction):
    candidates = state["candidates"]
    if candidates:
      index = next((i for i, c in enumerate(candidates) if c is state["selected"]), -1 if direction > 0 else 0)
      attempt(_global_fit_choose, state, candidates[(index+direction) % len(candidates)])
      update()
  def start(*_args):
    attempt(_global_fit_start, state, selection.get_active_id(), fields[0].get_text(), fields[1].get_text(),
            blur_b.get_text() if blur.get_active() else 0)
  def original(*_args):
    attempt(_global_fit_choose, state)
    update()
  state["update"] = update
  def on_result_key(_controller, keyval, _keycode, _modifiers):
    if keyval in (65362, 65364) and not state["running"] and rows.get_selected_row() is not None:
      navigate(-1 if keyval == 65362 else 1)
      row = rows.get_selected_row()
      if row is not None:
        row.grab_focus()
      return True
    return False
  keys = Gtk.EventControllerKey()
  keys.set_propagation_phase(Gtk.PropagationPhase.CAPTURE)
  keys.connect("key-pressed", on_result_key)
  rows.add_controller(keys)
  rows.connect("row-selected", choose)
  actions["start"].connect("clicked", start)
  actions["move"].connect("clicked", lambda *_: attempt(_global_fit_move_to_center, state))
  actions["reset"].connect("clicked", lambda *_: attempt(_global_fit_reset, state))
  actions["stop"].connect("clicked", lambda *_: attempt(_global_fit_stop, state))
  actions["previous"].connect("clicked", lambda *_: navigate(-1))
  actions["next"].connect("clicked", lambda *_: navigate(1))
  actions["original"].connect("clicked", original)
  actions["accept"].connect("clicked", lambda *_: attempt(_global_fit_finish, state, True))
  actions["close"].connect("clicked", lambda *_: attempt(_global_fit_finish, state))
  window.connect("close-request", lambda *_: attempt(_global_fit_finish, state) or True)
  window.set_child(box)
  state["widgets"] = dict(actions, rows=rows, scope=selection, trials=fields[0], radius=fields[1],
                          progress=progress, blur=blur, blur_b=blur_b, heading=heading, keys=keys, identifier=label)
  _GLOBAL_FIT_STATE = state
  update()
  status.set_text("Uses the active scroll-wheel map, captured when this panel opens. Coarsen it first if desired.")
  window.present()
  return state


#Jiggle-fits active chain to map
def jiggle_fit_active_chain():
  if (imol_refinement_map()==-1):
    info_dialog("You must set a refinement map!")
    return None
  chain_context = _active_chain_context_or_status()
  if not chain_context:
    return None
  fit_chain_to_map_by_random_jiggle(chain_context["mol_id"], chain_context["chain_id"], 1000, 0.1)

#Jiggle-fit active chain to B-smoothed map
def jiggle_fit_active_chain_smooth():
  if (imol_refinement_map()==-1):
    info_dialog("You must set a refinement map!")
    return None
  chain_context = _active_chain_context_or_status()
  if not chain_context:
    return None
  original_map = imol_refinement_map()
  existing = set(map_molecule_list())
  blur = getattr(coot, "sharpen_blur_map", None)
  if not callable(blur):
    info_dialog("This Coot build cannot create a disposable blurred map.")
    return None
  blurred = blur(original_map, 200)
  if blurred in existing or blurred not in map_molecule_list():
    info_dialog("Could not create a blurred fitting map; original map unchanged.")
    return None
  try:
    set_map_displayed(blurred, 0)
    set_imol_refinement_map(blurred)
    fit_chain_to_map_by_random_jiggle(chain_context["mol_id"], chain_context["chain_id"], 1000, 0.1)
  finally:
    set_imol_refinement_map(original_map)
    close_molecule(blurred)
    
#Jiggle-fits active chain to map (more thorough)
def jiggle_fit_active_chain_slow():
  if (imol_refinement_map()==-1):
    info_dialog("You must set a refinement map!")
    return None
  chain_context = _active_chain_context_or_status()
  if not chain_context:
    return None
  score_strings=[]
  for scale_fac in range(1,1000,100):
    scale_fac2=scale_fac/10.0
    score=fit_chain_to_map_by_random_jiggle(chain_context["mol_id"], chain_context["chain_id"], 100, scale_fac2)
    score_strings.append("scale:{scale} score:{score}".format(scale=scale_fac2, score=score))
  print(" ".join(score_strings))
    
#Fits all polymer chains to map
def jiggle_fit_all_chains():
  if (imol_refinement_map()==-1):
    info_dialog("You must set a refinement map!")
    return None
  mol_id = _active_molecule_or_status()
  if mol_id is None:
    return None
  with _grouped_model_edit(mol_id):
    for ch_id in chain_ids(mol_id):
      if is_polymer(mol_id,ch_id)==1:
        fit_chain_to_map_by_random_jiggle(mol_id,ch_id,1000,0.1)
        accept_regularizement()

#Jiggle-fit current molecule to map
def jiggle_fit_active_mol():
  if (imol_refinement_map()==-1):
    info_dialog("You must set a refinement map!")
    return None
  mol_id = _active_molecule_or_status()
  if mol_id is None:
    return None
  fit_molecule_to_map_by_random_jiggle(mol_id,1000,0.1)
    
#Clear labels and distances
def clear_distances_and_labels():
  remove_all_atom_labels()
  clear_measure_distances()
  
def _snap_click_to_nearby_segment_terminus(mol_id, ch_id, resno, max_distance=2):
  if mol_id not in model_molecule_list() or not ch_id or resno is False:
    return None
  try:
    first_in_seg = first_residue_in_seg(mol_id, ch_id, resno)
    last_in_seg = last_residue_in_seg(mol_id, ch_id, resno)
  except Exception:
    return None

  candidates = [
    ("N", first_in_seg, abs(first_in_seg - resno)),
    ("C", last_in_seg, abs(last_in_seg - resno)),
  ]
  candidates = [candidate for candidate in candidates if candidate[2] <= max_distance]
  if not candidates:
    return None
  candidates.sort(key=lambda candidate: candidate[2])
  terminus_kind, terminus_resno, _distance = candidates[0]
  return (terminus_kind, terminus_resno)


#click the start and end point, then fit the gap between them with polyala
def fit_polyala_gui():
  def fit_polyala(res1,res2):
    mol_id_1=_click_spec_imol(res1)
    mol_id_2=_click_spec_imol(res2)
    ch_id_1=_click_spec_chain_id(res1)
    ch_id_2=_click_spec_chain_id(res2)
    resno_1=_click_spec_res_no(res1)
    resno_2=_click_spec_res_no(res2)
    if (mol_id_1!=mol_id_2) or (ch_id_1!=ch_id_2):
      info_dialog("Start and end residues must be in the same molecule and chain!")
      return None
    terminus_1 = _snap_click_to_nearby_segment_terminus(mol_id_1, ch_id_1, resno_1)
    terminus_2 = _snap_click_to_nearby_segment_terminus(mol_id_1, ch_id_1, resno_2)
    if terminus_1 is None or terminus_2 is None:
      info_dialog("Click two termini flanking the loop, or within two residues of them.")
      return None

    clicked_termini = sorted([terminus_1, terminus_2], key=lambda item: item[1])
    lower_kind, lower_resno = clicked_termini[0]
    upper_kind, upper_resno = clicked_termini[1]
    if lower_kind != "C" or upper_kind != "N":
      info_dialog("Clicks must identify the C-terminus and N-terminus flanking the loop.")
      return None

    gap_start=lower_resno + 1
    gap_end=upper_resno - 1
    length=gap_end-gap_start+1
    if length <= 0:
      info_dialog("The selected termini do not flank a missing loop.")
      return None
    loop_seq=length*"A"
    fit_gap(mol_id_1,ch_id_1,gap_start,gap_end,loop_seq,1)
  _script_user_defined_click(2, fit_polyala,
    prompt="Build poly-Ala loop: click the two termini flanking the gap. Tab cancels.")
  
# Try to rebuild with db_mainchain after fit_gap?

#Real space refine for keyboard shortcut
def refine_click():
  def refine_click_a(res1,res2):
    mol_id_1=_click_spec_imol(res1)
    mol_id_2=_click_spec_imol(res2)
    ch_id_1=_click_spec_chain_id(res1)
    ch_id_2=_click_spec_chain_id(res2)
    resno_1=_click_spec_res_no(res1)
    resno_2=_click_spec_res_no(res2)
    alt_conf_1=_click_spec_alt_conf(res1)
    alt_conf_2=_click_spec_alt_conf(res2)
    if _clicked_same_chain_residue_range(res1, res2) is None:
      return
    if resno_2>=resno_1:
      refine_zone(mol_id_1,ch_id_1,resno_1,resno_2,alt_conf_1)
    else:
      refine_zone(mol_id_1,ch_id_1,resno_2,resno_1,alt_conf_1)
  _script_user_defined_click(2, refine_click_a,
    prompt="Refine range: click the start and end residues. Tab cancels.")
  
#Refine a clicked residue range.
def refine_residues_click():
  def refine_residues_click_a(res1,res2):
    mol_id_1=_click_spec_imol(res1)
    mol_id_2=_click_spec_imol(res2)
    ch_id_1=_click_spec_chain_id(res1)
    ch_id_2=_click_spec_chain_id(res2)
    resno_1=_click_spec_res_no(res1)
    resno_2=_click_spec_res_no(res2)
    alt_conf_1=_click_spec_alt_conf(res1)
    if imol_refinement_map()==-1:
      add_status_bar_text("You need to set a refinement map")
      return
    if (mol_id_1==mol_id_2) and (ch_id_1==ch_id_2):
      if resno_1>=resno_2:
        resno_1,resno_2=resno_2,resno_1 #Swap res numbers if wrong way round
      previous_immediate_replacement = refinement_immediate_replacement_state()
      try:
        set_refinement_immediate_replacement(0)
        refine_zone(mol_id_1,ch_id_1,resno_1,resno_2,alt_conf_1)
      finally:
        set_refinement_immediate_replacement(previous_immediate_replacement)
    else:
      info_dialog("Sorry, start and end residues must be in same mol and chain!")
  _script_user_defined_click(2, refine_residues_click_a,
    prompt="Refine range: click the start and end residues. Tab cancels.")
  
def refine_residues_range(mol_id,ch_id,resno_1,resno_2,ins_code_1):
  if resno_1>=resno_2:
    resno_1,resno_2=resno_2,resno_1 #Swap res numbers if wrong way round
  res_list=[]
  for resn in range(resno_1,resno_2+1):
    res_triple=[ch_id,resn,ins_code_1] #Does not account for varying ins code. Should be able to fix using residues_matching_criteria(), e.g.:
    #residues_matching_criteria(0, lambda chain_id,resno,ins_code,serial: True), except substitute a test func for true
    #that evaluates to true if resno, mol id and ch_id match, and returns ins_code (third item in output triple)
    res_list.append(res_triple) #append rather than adding, bc making list of lists
  refine_residues(mol_id,res_list)

  
def _clicked_same_chain_residue_range(res1, res2, allow_single_residue=True, error_message="Start and end residues must be in the same molecule and chain!"):
  mol_id_1=_click_spec_imol(res1)
  mol_id_2=_click_spec_imol(res2)
  ch_id_1=_click_spec_chain_id(res1)
  ch_id_2=_click_spec_chain_id(res2)
  resno_1=_click_spec_res_no(res1)
  resno_2=_click_spec_res_no(res2)
  if mol_id_1 < 0 or (mol_id_1!=mol_id_2) or (ch_id_1!=ch_id_2):
    info_dialog(error_message)
    return None
  if (resno_1, _click_spec_ins_code(res1)) == (resno_2, _click_spec_ins_code(res2)) and not allow_single_residue:
    info_dialog(error_message)
    return None
  return (mol_id_1, ch_id_1, min(resno_1, resno_2), max(resno_1, resno_2))


def _clicked_residue_specs(res1, res2):
  clicked = _clicked_same_chain_residue_range(res1, res2)
  if clicked is None:
    return []
  mol_id, chain, _first, _last = clicked
  specs = _tq_chain_residue_specs(mol_id, chain)
  first = (chain, _click_spec_res_no(res1), _click_spec_ins_code(res1))
  last = (chain, _click_spec_res_no(res2), _click_spec_ins_code(res2))
  if first not in specs or last not in specs:
    info_dialog("One of the clicked residues is no longer present.")
    return []
  a, b = sorted((specs.index(first), specs.index(last)))
  return specs[a:b+1]


def _atom_selection_for_residue_range(ch_id, res_start, res_end):
  if res_start == res_end:
    return "//%s/%s" % (ch_id, res_start)
  return "//%s/%s-%s" % (ch_id, res_start, res_end)


  
#Copy fragment (click start and end)
def copy_frag_by_click():
  def copy_frag(res1,res2):
    clicked_range = _clicked_same_chain_residue_range(res1, res2)
    if clicked_range is None:
      return None
    mol_id, ch_id, res_start, res_end = clicked_range
    specs = _clicked_residue_specs(res1, res2)
    if specs:
      return _copy_residue_specs_verified(mol_id, specs)
  _script_user_defined_click(2, copy_frag,
    prompt="Copy fragment: click the start and end residues. Tab cancels.")
  
#Cut fragment (click start and end)
def cut_frag_by_click():
  def cut_frag(res1,res2):
    clicked_range = _clicked_same_chain_residue_range(res1, res2)
    if clicked_range is None:
      return None
    mol_id, ch_id, res_start, res_end = clicked_range
    specs = _clicked_residue_specs(res1, res2)
    if specs:
      return _cut_residue_specs_verified(mol_id, specs)
  _script_user_defined_click(2, cut_frag,
    prompt="Cut fragment: click the start and end residues. Tab cancels.")


#Delete sidechain range (click start and end)
def delete_sidechain_range(mol_id, ch_id, res_start, res_end):
  for resno in range(res_start, res_end + 1):
    coot.delete_residue_sidechain(mol_id, ch_id, resno, "", 0)

# Sequence-register editing groups native backups and preserves rotamer settings.
_REGISTER_AA = dict(zip("ARNDCQEGHILKMFPSTWYV", (
  "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
  "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL")))
_REGISTER_PROTEIN_NAMES = set(_REGISTER_AA.values()) | {"MSE", "UNK"}


def _register_clean_sequence(text):
  lines = str(text).strip().splitlines()
  if sum(line.lstrip().startswith(">") for line in lines) > 1:
    raise ValueError("Supply one sequence, not multiple FASTA records.")
  sequence = "".join("".join(line.split()) for line in lines if not line.lstrip().startswith(">")).upper()
  if not sequence or set(sequence) - (set(_REGISTER_AA) | {"X"}):
    raise ValueError("Use a protein sequence or single FASTA record. Gaps, stop symbols and nonstandard letters are unsupported; X is allowed outside the target region.")
  return sequence


def _register_chain_snapshot(mol_id, chain_id):
  if not valid_model_molecule_qm(mol_id):
    raise ValueError("The target model is no longer available.")
  if coot.n_models(mol_id) != 1:
    raise ValueError("Register shifting currently requires a single-model molecule.")
  records = []
  for ch, number, ins in _all_residue_specs_for_colouring(mol_id):
    if ch != chain_id:
      continue
    atoms = coot.residue_info_py(mol_id, ch, number, ins) or []
    signature = tuple(sorted((a[0][0], a[0][1], tuple(float(x) for x in a[2]), repr(a[1])) for a in atoms))
    alternatives = sorted({a[0][1] for a in atoms if a[0][1]})
    preferred = max(alternatives, key=lambda alt: sum(float(a[1][0]) for a in atoms if a[0][1] == alt) /
                    sum(a[0][1] == alt for a in atoms)) if alternatives else ""
    shared = {a[0][0] for a in atoms if not a[0][1]}
    backbone = {a[0][0].strip(): tuple(float(x) for x in a[2]) for a in atoms
                if a[0][0].strip() in {"N", "CA", "C", "O"} and
                (not a[0][1] or (a[0][1] == preferred and a[0][0] not in shared))}
    records.append({"number": number, "ins": ins, "name": residue_name(mol_id, ch, number, ins).strip(),
                    "alt": any(bool(a[0][1]) for a in atoms), "backbone": backbone, "atoms": signature})
  return sorted(records, key=lambda r: (r["number"], r["ins"]))


def _register_connected(left, right):
  if (left["name"] not in _REGISTER_PROTEIN_NAMES or right["name"] not in _REGISTER_PROTEIN_NAMES
      or left["ins"] or right["ins"] or right["number"] != left["number"] + 1):
    return False
  carbon, nitrogen = left["backbone"].get("C"), right["backbone"].get("N")
  return carbon is not None and nitrogen is not None and sum((a-b)**2 for a, b in zip(carbon, nitrogen)) <= 2.2**2


def _register_plan(records, start, end, offset, sequence, sequence_start):
  """Pure preflight: trim only the nearest neighbour's overlapping terminus."""
  start, end = sorted((int(start), int(end)))
  offset, sequence_start = int(offset), int(sequence_start)
  if not offset:
    raise ValueError("Enter a nonzero register shift; +2 means residue 100 becomes 102.")
  sequence = _register_clean_sequence(sequence)
  selected = [r for r in records if start <= r["number"] <= end]
  if [r["number"] for r in selected] != list(range(start, end+1)):
    raise ValueError("The selected range must have unique consecutive residue numbers, with no missing residues or insertion codes.")
  for r in selected:
    if r["ins"] or r["name"] not in _REGISTER_PROTEIN_NAMES or set(r["backbone"]) != {"N", "CA", "C", "O"}:
      raise ValueError("The selected range requires protein residues with complete N/CA/C/O atoms and no insertion codes.")
  if any(not _register_connected(a, b) for a, b in zip(selected, selected[1:])):
    raise ValueError("The selected range crosses a backbone break. Select a single connected region.")
  first, last = start+offset, end+offset
  if first < sequence_start or last >= sequence_start+len(sequence):
    raise ValueError("The shifted region is outside the supplied sequence numbering. Check the sequence's first residue number.")
  target = sequence[first-sequence_start:last-sequence_start+1]
  if set(target) - set(_REGISTER_AA):
    raise ValueError("The target region contains X; specify those residue identities before shifting.")
  outsiders = [r for r in records if not start <= r["number"] <= end]
  trim = [r for r in outsiders if first <= r["number"] <= last]
  if trim:
    side = [r for r in records if r["number"] > end] if offset > 0 else [r for r in records if r["number"] < start]
    if offset < 0:
      side = list(reversed(side))
    neighbour = []
    for r in side:
      if neighbour:
        a, b = (neighbour[-1], r) if offset > 0 else (r, neighbour[-1])
        if not _register_connected(a, b):
          break
      neighbour.append(r)
    trimmed_numbers = {(r["number"], r["ins"]) for r in trim}
    prefix = neighbour[:len(trim)]
    if trimmed_numbers != {(r["number"], r["ins"]) for r in prefix}:
      raise ValueError("The shift collides with more than the adjacent region's terminus.")
    if len(trim) >= len(neighbour):
      raise ValueError("The shift would delete the entire neighbouring region. Use a smaller shift or handle that region separately.")
    if any(r["ins"] or r["name"] not in _REGISTER_PROTEIN_NAMES for r in trim):
      raise ValueError("Overlaps involving ligands or insertion codes cannot be trimmed automatically.")
  remaining = [r for r in outsiders if r not in trim]
  if any((r["number"] < start and r["number"] >= first) or
         (r["number"] > end and r["number"] <= last) for r in remaining):
    raise ValueError("The shift would move the region past another residue or fragment.")
  changes = []
  for r, letter in zip(selected, target):
    new_name = _REGISTER_AA[letter]
    equivalent = r["name"] == new_name or (r["name"] == "MSE" and new_name == "MET")
    changes.append((r["number"], r["number"]+offset, r["name"], r["name"] if equivalent else new_name, not equivalent or r["alt"]))
  return {"start": start, "end": end, "offset": offset, "target": target,
          "sequence": sequence, "sequence_start": sequence_start, "trim": trim, "changes": changes,
          "snapshot": records, "selected": selected}


def _register_preview_text(plan, map_id):
  reverse = {value: key for key, value in _REGISTER_AA.items()}
  old = "".join(reverse.get(r["name"], "M" if r["name"] == "MSE" else "X") for r in plan["selected"])
  first, last = plan["start"]+plan["offset"], plan["end"]+plan["offset"]
  changed = sum(c[4] for c in plan["changes"])
  text = [f"{plan['start']}-{plan['end']}  ->  {first}-{last}    ({plan['offset']:+d})", "",
          f"Before  {plan['start']}-{old}-{plan['end']}",
          f"After   {first}-{plan['target']}-{last}", "",
          f"Rebuild/fit: {changed} of {len(plan['selected'])} residues | Map #{map_id}",
          "Delete: " + (", ".join(f"{r['number']} {r['name']}" for r in plan["trim"]) or "none"),
          "", "Current rotamer/backrub settings apply.",
          "Vacated positions stay gaps. One Undo restores the whole shift."]
  collapsed = [str(r["number"]) for r in plan["selected"] if r["alt"]]
  if collapsed:
    text.insert(7, "Altlocs -> single rotamer: " + ", ".join(collapsed))
  return "\n".join(text)


def _register_collapse_alternates(mol_id, chain_id, number):
  atoms = coot.residue_info_py(mol_id, chain_id, number, "") or []
  alternatives = sorted({a[0][1] for a in atoms if a[0][1]})
  if not alternatives:
    return
  preferred = max(alternatives, key=lambda alt: sum(float(a[1][0]) for a in atoms if a[0][1] == alt) /
                  sum(a[0][1] == alt for a in atoms))
  shared = {a[0][0] for a in atoms if not a[0][1]}
  occupancy = {}
  for atom in atoms:
    occupancy[atom[0][0]] = occupancy.get(atom[0][0], 0.0) + float(atom[1][0])
  for atom in atoms:
    name, alt = atom[0]
    if alt and (alt != preferred or name in shared):
      coot.delete_atom(mol_id, chain_id, number, "", name, alt)
  for atom in atoms:
    name, alt = atom[0]
    if alt == preferred and name not in shared:
      # Native implementations can return 0 even after successfully editing.
      coot.set_atom_string_attribute(mol_id, chain_id, number, "", name, alt, "alt-conf", "")
      updated = coot.residue_info_py(mol_id, chain_id, number, "") or []
      if not any(a[0][0] == name and not a[0][1] for a in updated) or any(
          a[0][0] == name and a[0][1] == alt for a in updated):
        raise RuntimeError(f"Could not clear alternate identifier for {number}:{name}.")
      coot.set_atom_attribute(mol_id, chain_id, number, "", name, "", "occ", min(1.0, occupancy[name]))
  remaining = coot.residue_info_py(mol_id, chain_id, number, "") or []
  if any(a[0][1] for a in remaining) or len({a[0][0] for a in remaining}) != len(remaining):
    raise RuntimeError(f"Could not safely collapse alternate conformers for residue {number}.")


def _apply_register_shift(mol_id, chain_id, plan, map_id, anchor_number=None):
  global NAVIGATION_LAST_RESIDUE
  try:
    if _register_chain_snapshot(mol_id, chain_id) != plan["snapshot"]:
      raise ValueError("The chain changed after preview. Preview again before applying.")
    if map_id not in map_molecule_list():
      raise ValueError("The selected refinement map is no longer available.")
    if coot.backup_state(mol_id) != 1:
      raise ValueError("Enable Coot backups for this molecule before shifting its register.")
  except (ValueError, TypeError) as error:
    info_dialog(str(error))
    return False
  failed_fits = []
  stage = "backup preparation"
  try:
    # One native history entry covers trimming, renumbering and sidechain edits.
    coot.make_backup(mol_id)
    coot.turn_off_backup(mol_id)
    stage = "trimming"
    for r in plan["trim"]:
      coot.delete_residue(mol_id, chain_id, r["number"], r["ins"])
      if coot.residue_info_py(mol_id, chain_id, r["number"], r["ins"]):
        raise RuntimeError(f"Coot could not delete overlapping residue {r['number']}.")
    stage = "renumbering"
    if not coot.renumber_residue_range(mol_id, chain_id, plan["start"], plan["end"], plan["offset"]):
      raise RuntimeError("Coot rejected renumbering.")
    stage = "mutation/autofitting"
    for old_no, new_no, old_name, new_name, changed in plan["changes"]:
      if changed:
        _register_collapse_alternates(mol_id, chain_id, new_no)
        if not coot.mutate(mol_id, chain_id, new_no, "", new_name):
          raise RuntimeError(f"Coot could not mutate residue {new_no} to {new_name}.")
        if new_name not in {"ALA", "GLY"}:
          score = coot.auto_fit_best_rotamer(mol_id, chain_id, new_no, "", "", map_id, 1, 0.3)
          if not isinstance(score, (int, float)) or not math.isfinite(score) or score <= -999:
            failed_fits.append(new_no)
  except Exception as error:
    NAVIGATION_LAST_RESIDUE = None
    info_dialog(f"Register shift stopped during {stage}. Earlier edits remain; use Coot Undo to recover.\n\n{error}")
    return False
  finally:
    coot.turn_on_backup(mol_id)
    _mark_model_changed(mol_id)
  NAVIGATION_LAST_RESIDUE = None
  destination = (anchor_number if anchor_number is not None and plan["start"] <= anchor_number <= plan["end"] else plan["start"]) + plan["offset"]
  coot.set_go_to_atom_molecule(mol_id)
  coot.set_go_to_atom_chain_residue_atom_name_full(chain_id, destination, "", " CA ", "")
  message = f"Shifted {chain_id}:{plan['start']}-{plan['end']} by {plan['offset']:+d}; trimmed {len(plan['trim'])}, mutated {sum(c[4] for c in plan['changes'])}."
  print(message)
  add_status_bar_text(_navigation_status_bar_label(mol_id, chain_id, destination, "") + " | " + message)
  if failed_fits:
    info_dialog(message + "\n\nRotamer fitting could not complete for: " + ", ".join(map(str, failed_fits)))
  return True


def _register_sequence_origin(records, sequence):
  """Infer a numbering offset only from a strong, unique identity match."""
  sequence = _register_clean_sequence(sequence)
  reverse = {value: key for key, value in _REGISTER_AA.items()}
  reverse["MSE"] = "M"
  known = [(r["number"], reverse[r["name"]]) for r in records
           if not r["ins"] and r["name"] in reverse]
  if len(known) < 3:
    return None
  candidates = set()
  for triple in zip(known, known[1:], known[2:]):
    if triple[2][0] != triple[0][0] + 2 or triple[1][0] != triple[0][0] + 1:
      continue
    motif = "".join(letter for _, letter in triple)
    position = sequence.find(motif)
    while position >= 0:
      candidates.add(triple[0][0] - position)
      position = sequence.find(motif, position + 1)
  scores = [(sum(0 <= number-origin < len(sequence) and sequence[number-origin] == letter
                 for number, letter in known), origin) for origin in candidates]
  scores.sort(reverse=True)
  if not scores or scores[0][0] < max(3, 0.8 * len(known)):
    return None
  if len(scores) > 1 and scores[0][0] == scores[1][0]:
    return None
  return scores[0][1]


def _register_reference_sequence(mol_id, chain_id, records):
  """Prefer reference sequence; coordinates are a last-resort, numbered fallback."""
  getter = getattr(coot, "sequence_info_py", None)
  if callable(getter):
    associated = [item[1] for item in getter(mol_id) or [] if len(item) == 2 and item[0] == chain_id]
    if len(associated) == 1 and associated[0]:
      return _register_clean_sequence(associated[0]), "Chain-associated reference sequence loaded."
  path = _annotation_source_path(mol_id)
  if path:
    names = []
    with open(path, encoding="utf-8", errors="replace") as handle:
      for line in handle:
        if line.startswith("SEQRES") and line[11:12].strip() == chain_id:
          names.extend(line[19:70].split())
    if names:
      reverse = {value: key for key, value in _REGISTER_AA.items()}
      reverse["MSE"] = "M"
      return "".join(reverse.get(name, "X") for name in names), "PDB SEQRES reference sequence loaded."
  reverse = {value: key for key, value in _REGISTER_AA.items()}
  reverse["MSE"] = "M"
  numbered = {r["number"]: reverse[r["name"]] for r in records
              if not r["ins"] and r["name"] in reverse}
  if numbered and max(numbered)-min(numbered) < 100000:
    return "".join(numbered.get(i, "X") for i in range(min(numbered), max(numbered)+1)), (
      "Model-derived sequence loaded (not an independent reference). Missing positions are X; paste a reference if existing identities are wrong.")
  return "", "Paste the reference sequence or a single FASTA record."


def _register_shift_model_signature(mol_id):
  """Guard native Undo, allowing the precision of Coot's PDB backups."""
  digest = hashlib.sha256()
  def rounded(value, digits):
    value = round(float(value), digits)
    return value if value else 0.0
  for spec, record in sorted(_model_coordinate_snapshot(mol_id).items()):
    atoms = sorted((a[0][0], a[0][1], rounded(a[1][0], 2), rounded(a[1][1], 2),
                    str(a[1][2]).strip(), tuple(rounded(x, 3) for x in a[2]))
                   for a in record["atoms"])
    digest.update(repr((spec, record["name"], atoms)).encode("utf-8"))
  return digest.digest()


def _register_shift_dialog(mol_id, chain_id, start, end, anchor_number=None, selection_cleanup=None):
  required = ("n_models", "residue_info_py", "renumber_residue_range", "delete_residue", "mutate",
              "auto_fit_best_rotamer", "backup_state", "make_backup", "turn_off_backup", "turn_on_backup",
              "delete_atom", "set_atom_string_attribute",
              "set_atom_attribute", "set_go_to_atom_molecule",
              "set_go_to_atom_chain_residue_atom_name_full")
  missing = [name for name in required if not callable(getattr(coot, name, None))]
  if missing or globals().get("Gtk") is None:
    if selection_cleanup:
      selection_cleanup()
    info_dialog("Register shift is unavailable in this build. " + ", ".join(missing))
    return
  records = _register_chain_snapshot(mol_id, chain_id)
  default_sequence, sequence_message = _register_reference_sequence(mol_id, chain_id, records)
  origin = None
  if default_sequence:
    try:
      origin = _register_sequence_origin(records, default_sequence)
      if sequence_message.startswith("Model-derived"):
        origin = min(r["number"] for r in records if not r["ins"] and r["name"] in _REGISTER_PROTEIN_NAMES - {"UNK"})
    except (ValueError, TypeError):
      pass
  window = Gtk.Window()
  def close_dialog(*_args):
    if selection_cleanup:
      selection_cleanup()
    window.destroy()
  def closing(*_args):
    if selection_cleanup:
      selection_cleanup()
    return False
  window.connect("close-request", closing)
  region = {"start": min(start, end), "end": max(start, end), "anchor": anchor_number}
  history = []
  def update_title():
    window.set_title(f"Shift register: molecule #{mol_id}, {chain_id}:{region['start']}-{region['end']}")
  update_title()
  window.set_default_size(620, 440)
  box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=8)
  for side in ("start", "end", "top", "bottom"):
    getattr(box, "set_margin_" + side)(12)
  source = ("PDB SEQRES" if sequence_message.startswith("PDB") else
            "Model-derived (verify identities)" if sequence_message.startswith("Model-derived") else
            "Associated sequence" if default_sequence else "Reference sequence required")
  needs_sequence_input = not default_sequence or origin is None or sequence_message.startswith("Model-derived")
  note = Gtk.Label(label=source)
  note.set_tooltip_text(sequence_message + "\nPositive shift: 100 + 2 becomes 102.")
  note.set_wrap(True)
  note.set_xalign(0)
  box.append(note)
  entries = []
  origin_row = None
  for label, value in (("Shift (+/-)", "1"), ("Sequence starts at", str(origin) if origin is not None else "")):
    row = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
    row.append(Gtk.Label(label=label))
    entry = Gtk.Entry()
    entry.set_text(value)
    row.append(entry)
    entries.append(entry)
    if label == "Sequence starts at":
      origin_row = row
    else:
      box.append(row)
  sequence_view = Gtk.TextView()
  sequence_view.set_wrap_mode(Gtk.WrapMode.WORD_CHAR)
  sequence_view.get_buffer().set_text(default_sequence)
  sequence_scroll = Gtk.ScrolledWindow()
  sequence_scroll.set_min_content_height(110)
  sequence_scroll.set_child(sequence_view)
  reference = Gtk.Expander(label="Reference sequence needs verification")
  reference.set_visible(needs_sequence_input)
  reference.set_expanded(needs_sequence_input)
  sequence_controls = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6)
  sequence_controls.append(origin_row)
  sequence_controls.append(sequence_scroll)
  reference.set_child(sequence_controls)
  box.append(reference)
  detail = Gtk.TextView()
  detail.set_editable(False)
  detail.set_monospace(True)
  detail.set_wrap_mode(Gtk.WrapMode.NONE)
  detail_scroll = Gtk.ScrolledWindow()
  detail_scroll.set_vexpand(True)
  detail_scroll.set_child(detail)
  box.append(detail_scroll)
  buttons = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
  apply_button = Gtk.Button(label="Apply register shift")
  apply_button.set_sensitive(False)
  undo_button = Gtk.Button(label="Undo")
  undo_button.set_sensitive(False)
  preview_button = Gtk.Button(label="Preview")
  undo_available = all(callable(getattr(coot, name, None)) for name in ("set_undo_molecule", "apply_undo"))
  undo_button.set_tooltip_text("Undo the last shift in this dialog, including trimmed residues."
                               if undo_available else "Native Undo is unavailable in this Coot build.")
  state = {}
  blocked = False
  def invalidate(*_args):
    state.clear()
    apply_button.set_sensitive(False)
    detail.get_buffer().set_text("Preview required after changing sequence or numbering.")
  for entry in entries:
    entry.connect("changed", invalidate)
  sequence_view.get_buffer().connect("changed", invalidate)
  def preview(*_args):
    invalidate()
    try:
      if blocked:
        raise ValueError("Reopen the shift dialog before making further edits.")
      if history and not history[-1]["success"]:
        raise ValueError("Undo the incomplete shift before trying another register.")
      pending = getattr(coot, "get_continue_updating_refinement_atoms_state", None)
      if callable(pending) and pending():
        raise ValueError("Accept or cancel the current refinement before shifting registers.")
      buffer = sequence_view.get_buffer()
      sequence = buffer.get_text(buffer.get_start_iter(), buffer.get_end_iter(), False)
      snapshot = _register_chain_snapshot(mol_id, chain_id)
      offset = int(entries[0].get_text())
      try:
        clean = _register_clean_sequence(sequence)
        sequence_start = int(entries[1].get_text())
        first, last = region["start"]+offset, region["end"]+offset
        if first < sequence_start or last >= sequence_start+len(clean):
          raise ValueError("The shifted range is outside the reference sequence. Check the sequence and its start.")
        if set(clean[first-sequence_start:last-sequence_start+1]) - set(_REGISTER_AA):
          raise ValueError("The target sequence has unknown residues; enter their identities.")
      except (ValueError, TypeError):
        reference.set_visible(True)
        reference.set_expanded(True)
        raise
      plan = _register_plan(snapshot, region["start"], region["end"], offset, clean, sequence_start)
      map_id = imol_refinement_map()
      if map_id not in map_molecule_list():
        raise ValueError("Select a refinement map in Coot before previewing.")
      state.update(plan=plan, map_id=map_id)
      detail.get_buffer().set_text(_register_preview_text(plan, map_id))
      apply_button.set_sensitive(True)
    except (ValueError, TypeError) as error:
      detail.get_buffer().set_text(str(error))
  def apply(*_args):
    if not state:
      return
    plan, map_id = state["plan"], state["map_id"]
    invalidate()
    try:
      before = _register_shift_model_signature(mol_id)
      try:
        success = _apply_register_shift(mol_id, chain_id, plan, map_id, region["anchor"])
      except Exception as error:
        traceback.print_exc()
        success = False
        info_dialog("Register shift stopped; Undo any partial edits.\n\n" + str(error))
      after = _register_shift_model_signature(mol_id)
      if after != before:
        history.append({"before": before, "after": after, "region": dict(region), "success": success})
      undo_button.set_sensitive(bool(history) and undo_available)
      if success:
        region["start"] += plan["offset"]
        region["end"] += plan["offset"]
        if region["anchor"] is not None:
          region["anchor"] += plan["offset"]
        update_title()
        if selection_cleanup:
          selection_cleanup()
        detail.get_buffer().set_text("Applied. Undo to try from the previous register, or preview another shift.\n\n" +
                                     _register_preview_text(plan, map_id))
      else:
        detail.get_buffer().set_text("Shift not completed. Undo any partial edits before trying again.")
    except Exception as error:
      traceback.print_exc()
      detail.get_buffer().set_text(str(error))

  def undo(*_args):
    global NAVIGATION_LAST_RESIDUE
    nonlocal blocked
    if not history or not undo_available:
      return
    invalidate()
    try:
      if not valid_model_molecule_qm(mol_id):
        raise ValueError("The model has been closed.")
      displayed = getattr(coot, "mol_is_displayed", None)
      if callable(displayed) and not displayed(mol_id):
        raise ValueError("Display this model before undoing its register shift.")
      modified = getattr(coot, "have_unsaved_changes_p", None)
      if callable(modified) and modified(mol_id) != 1:
        raise ValueError("This model has no pending edits to undo.")
      pending = getattr(coot, "get_continue_updating_refinement_atoms_state", None)
      if callable(pending) and pending():
        raise ValueError("Accept or cancel the current refinement before undoing.")
      if _register_shift_model_signature(mol_id) != history[-1]["after"]:
        raise ValueError("The model changed outside this dialog. Use Coot's Undo for those edits first.")
      coot.set_undo_molecule(mol_id)
      if coot.apply_undo() != 1:
        raise ValueError("Coot could not undo this shift; check that the model is displayed and its backups are available.")
      _mark_model_changed(mol_id)
      NAVIGATION_LAST_RESIDUE = None
      restored = _register_shift_model_signature(mol_id)
      if restored != history[-1]["before"]:
        history.clear()
        blocked = True
        undo_button.set_sensitive(False)
        raise ValueError("Undo restored a different model state. Reopen the shift dialog before making further edits.")
      region.update(history.pop()["region"])
      update_title()
      undo_button.set_sensitive(bool(history))
      if selection_cleanup:
        selection_cleanup()
      destination = region["anchor"] if region["anchor"] is not None else region["start"]
      coot.set_go_to_atom_molecule(mol_id)
      coot.set_go_to_atom_chain_residue_atom_name_full(chain_id, destination, "", " CA ", "")
      add_status_bar_text(f"Undid register shift | {chain_id}:{region['start']}-{region['end']}")
      preview()
    except Exception as error:
      traceback.print_exc()
      detail.get_buffer().set_text(str(error))

  for label, action in (("Close", close_dialog), ("Preview", preview), ("Undo", undo)):
    if label == "Undo":
      button = undo_button
    elif label == "Preview":
      button = preview_button
    else:
      button = Gtk.Button(label=label)
    button.connect("clicked", action)
    buttons.append(button)
  apply_button.connect("clicked", apply)
  buttons.append(apply_button)
  box.append(buttons)
  window.set_child(box)
  window.present()
  window.register_shift_widgets = {"shift": entries[0], "origin": entries[1], "sequence": sequence_view,
                                   "reference": reference, "detail": detail, "apply": apply_button,
                                   "undo": undo_button, "preview": preview_button, "region": region}
  preview()
  return window


def shift_register_active_segment():
  active = _active_residue_or_status()
  if not active:
    return
  mol_id, chain_id, number, ins = active[:4]
  try:
    records = _register_chain_snapshot(mol_id, chain_id)
    index = next(i for i, r in enumerate(records) if r["number"] == number and r["ins"] == ins)
    left = right = index
    while left > 0 and _register_connected(records[left-1], records[left]):
      left -= 1
    while right+1 < len(records) and _register_connected(records[right], records[right+1]):
      right += 1
    _register_shift_dialog(mol_id, chain_id, records[left]["number"], records[right]["number"], number)
  except (ValueError, StopIteration) as error:
    info_dialog(f"Cannot select the active segment.\n\n{error}")


_REGISTER_SELECTION_MARKERS = []
_REGISTER_SELECTION_TOKEN = None


def _clear_register_selection_markers():
  had_markers = bool(_REGISTER_SELECTION_MARKERS)
  while _REGISTER_SELECTION_MARKERS:
    obj = _REGISTER_SELECTION_MARKERS.pop()
    try:
      coot.close_generic_object(obj)
    except Exception as error:
      print("Could not clear register endpoint marker:", error)
  draw = getattr(coot, "graphics_draw", None)
  if had_markers and callable(draw):
    draw()


def _mark_register_endpoint(spec, colour, markers=None, label="Register range endpoint"):
  markers = _REGISTER_SELECTION_MARKERS if markers is None else markers
  required = ("new_generic_object_number", "to_generic_object_add_point",
              "set_display_generic_object", "close_generic_object")
  if any(not callable(getattr(coot, name, None)) for name in required):
    return
  atoms = coot.residue_info_py(_click_spec_imol(spec), _click_spec_chain_id(spec),
                               _click_spec_res_no(spec), _click_spec_ins_code(spec)) or []
  candidates = []
  for name in ("CA", "C4'", "C4*", "P"):
    candidates = [a for a in atoms if a[0][0].strip() == name]
    if candidates:
      break
  if not candidates:
    return
  atom = max(candidates, key=lambda a: (not bool(a[0][1]), float(a[1][0])))
  obj = coot.new_generic_object_number(label)
  markers.append(obj)
  try:
    coot.to_generic_object_add_point(obj, colour, 8, *atom[2])
    coot.set_display_generic_object(obj, 1)
    draw = getattr(coot, "graphics_draw", None)
    if callable(draw):
      draw()
  except Exception as error:
    print("Could not draw register endpoint marker:", error)
    coot.close_generic_object(obj)
    markers.remove(obj)


def shift_register_clicked_range():
  global _REGISTER_SELECTION_TOKEN
  _clear_register_selection_markers()
  token = object()
  _REGISTER_SELECTION_TOKEN = token
  def cleanup():
    if _REGISTER_SELECTION_TOKEN is token:
      _clear_register_selection_markers()
  def show_click(spec, prompt):
    add_status_bar_text(_navigation_status_bar_label(
      _click_spec_imol(spec), _click_spec_chain_id(spec),
      _click_spec_res_no(spec), _click_spec_ins_code(spec)).replace("Current residue:", "Clicked residue:", 1) + " | " + prompt)
  def clicked(first, last):
    if _click_spec_ins_code(first) or _click_spec_ins_code(last):
      cleanup()
      info_dialog("Insertion-code endpoints are unsupported for register shifting.")
      return
    region = _clicked_same_chain_residue_range(first, last)
    if region:
      try:
        _register_shift_dialog(*region, anchor_number=_click_spec_res_no(first),
                               selection_cleanup=cleanup)
      except Exception:
        cleanup()
        raise
    else:
      cleanup()
  def show_endpoint(spec, index):
    show_click(spec, "Click the end residue (orange). Tab cancels." if index == 1 else "Range selected")
    _mark_register_endpoint(spec, "cyan" if index == 1 else "orange")
  _script_user_defined_click(2, clicked, on_pick=show_endpoint, cleanup=cleanup,
    prompt="Shift register: click start (cyan), then end (orange). Tab cancels.")


_REGISTER_REFERENCE_ORIGINS = {}


def _register_associate_reference(fragment, chain, sequence, origin):
  sequence = _register_clean_sequence(sequence)
  coot.delete_sequence_by_chain_id(fragment, chain)
  coot.assign_sequence_from_string(fragment, chain, sequence)
  associated = coot.sequence_info_py(fragment) or []
  if not any(item[0] == chain and _register_clean_sequence(item[1]) == sequence
             for item in associated if isinstance(item, (list, tuple)) and len(item) == 2):
    raise RuntimeError("Coot did not associate the reference sequence with the fragment.")
  _REGISTER_REFERENCE_ORIGINS[(fragment, chain)] = (sequence, int(origin))


def _register_terminal_residue_type(mol_id, chain, destination):
  cached = _REGISTER_REFERENCE_ORIGINS.get((mol_id, chain))
  if cached is not None:
    sequence, origin = cached
    associated = coot.sequence_info_py(mol_id) or []
    matches = [item[1] for item in associated if isinstance(item, (list, tuple))
               and len(item) == 2 and item[0] == chain]
    if len(matches) == 1 and _register_clean_sequence(matches[0]) == sequence:
      index = destination-origin
      if 0 <= index < len(sequence) and sequence[index] in _REGISTER_AA:
        return _REGISTER_AA[sequence[index]]
      add_status_bar_text("Reference sequence has no known residue at this position; using Coot's automatic type.")
    else:
      _REGISTER_REFERENCE_ORIGINS.pop((mol_id, chain), None)
  lookup = getattr(coot, "find_terminal_residue_type_py", None)
  if callable(lookup):
    residue_type = lookup(mol_id, chain, destination)
    if isinstance(residue_type, str) and residue_type in set(_REGISTER_AA.values()):
      return residue_type
  return "auto"


def _register_comparison_targets(records, centre, sequence, origin, shift_window=5, fragment_size=9):
  sequence = _register_clean_sequence(sequence)
  origin = int(origin)
  shift_window = int(shift_window)
  if shift_window < 0:
    raise ValueError("Shift window must be a non-negative integer.")
  fragment_size = int(fragment_size)
  if fragment_size < 1:
    raise ValueError("Fragment size must be a positive integer.")
  # Even lengths place the extra residue on the C-terminal side.
  start = centre-(fragment_size-1)//2
  end = start+fragment_size-1
  selected = [r for r in records if start <= r["number"] <= end]
  if [r["number"] for r in selected] != list(range(start, end+1)):
    raise ValueError(f"The fragment needs consecutive residues {start}-{end}, without insertion codes.")
  if any(r["ins"] or r["name"] not in _REGISTER_PROTEIN_NAMES or
         set(r["backbone"]) != {"N", "CA", "C", "O"} for r in selected):
    raise ValueError("The fragment must have complete protein backbone atoms and no insertion codes.")
  if any(not _register_connected(a, b) for a, b in zip(selected, selected[1:])):
    raise ValueError("The fragment crosses a backbone break.")
  if start-shift_window < origin or end+shift_window >= origin+len(sequence):
    raise ValueError(f"The reference sequence must cover residues {start-shift_window}-{end+shift_window}.")
  targets = []
  for offset in range(-shift_window, shift_window+1):
    target = sequence[start+offset-origin:end+offset-origin+1]
    if set(target) - set(_REGISTER_AA):
      raise ValueError("All candidate identities must be known; replace X in the reference sequence.")
    targets.append((offset, target))
  return start, end, targets


def _register_fit_fixed_backbone_rotamer(fragment, chain, number, map_id):
  """Score library conformers without the residue translations used by autofit."""
  rows = coot.score_rotamers_py(fragment, chain, number, "", "", map_id, 1, 0.3) or []
  candidates = []
  clash_acceptable = []
  for row in rows:
    if not isinstance(row, (list, tuple)) or len(row) < 5:
      raise RuntimeError("Unexpected rotamer scoring result.")
    density, clash = float(row[2]), float(row[4])
    if math.isfinite(density) and math.isfinite(clash) and clash >= 0:
      candidate = (density, row[0])
      candidates.append(candidate)
      if clash < 20:
        clash_acceptable.append(candidate)
  if not candidates:
    raise RuntimeError(f"No scorable fixed-backbone rotamer for {chain}:{number}.")
  # Unfitted neighbours can clash with every conformer; do not reject a register
  # before joint refinement has had a chance to resolve those contacts.
  best = max(clash_acceptable or candidates, key=lambda candidate: candidate[0])
  if not clash_acceptable:
    print(f"Register comparison {chain}:{number}: all rotamers clash; using best density fit.")
  def backbone_coordinates():
    return {tuple(atom[0]): tuple(float(value) for value in atom[2])
            for atom in coot.residue_info_py(fragment, chain, number, "") or []
            if atom[0][0].strip() in {"N", "CA", "C", "O"}}
  before = backbone_coordinates()
  if len(before) != 4:
    raise RuntimeError(f"Incomplete backbone at {chain}:{number}.")
  if not coot.set_residue_to_rotamer_name(fragment, chain, number, "", "", best[1]):
    raise RuntimeError(f"Could not set rotamer at {chain}:{number}.")
  after = backbone_coordinates()
  if before.keys() != after.keys() or any(
      not all(math.isfinite(value) for value in after[key]) or
      _distance_sq(before[key], after[key]) > 1e-10 for key in before):
    raise RuntimeError(f"Rotamer selection moved backbone atoms at {chain}:{number}; candidate excluded.")


def _register_sidechain_correlation(fragment, specs, target, map_id):
  """Mean residue CC beyond C-beta; Gly/Ala have no informative atoms here."""
  eligible = {tuple(spec) for spec, letter in zip(specs, target) if letter not in {"G", "A"}}
  scores = {}
  if eligible:
    # Mask 3 excludes main-chain atoms and CB in both CCP4 and Homebrew Coot.
    results = coot.map_to_model_correlation_per_residue_py(fragment, specs, 3, map_id) or []
    for item in results:
      if not isinstance(item, (list, tuple)) or len(item) < 2:
        continue
      spec = _coot_residue_spec_from_spec(item[0])
      if spec is None or tuple(spec) not in eligible:
        continue
      value = float(item[1])
      if math.isfinite(value) and -1 <= value <= 1:
        scores[tuple(spec)] = value
  return {"mean": sum(scores.values())/len(scores) if scores else None,
          "count": len(scores), "eligible": len(eligible)}


def _build_register_comparison_fragment(source, chain, start, end, offset, target, map_id, rigid_fit=False):
  """Build on a disposable copy; restore global refinement settings on every path."""
  fragment = coot.new_molecule_by_atom_selection(source, _atom_selection_for_residue_range(chain, start, end))
  if fragment == source or not valid_model_molecule_qm(fragment):
    raise RuntimeError("Coot could not create an independent fragment.")
  previous_map = coot.imol_refinement_map()
  previous_replacement = coot.refinement_immediate_replacement_state()
  try:
    coot.set_imol_refinement_map(map_id)
    # Optional rotamer post-refinement must auto-accept on the disposable copy.
    coot.set_refinement_immediate_replacement(1)
    coot.set_mol_displayed(fragment, 0)
    coot.set_molecule_name(fragment, f"Register {offset:+d}: {chain}:{start+offset}-{end+offset}")
    specs = [[chain, number+offset, ""] for number in range(start, end+1)]
    def fit_rotamers():
      for spec, letter in zip(specs, target):
        if _REGISTER_AA[letter] not in {"ALA", "GLY"}:
          _register_fit_fixed_backbone_rotamer(fragment, chain, spec[1], map_id)

    with _grouped_model_edit(fragment):
      if offset and not coot.renumber_residue_range(fragment, chain, start, end, offset):
        raise RuntimeError("Fragment renumbering failed.")
      for spec, letter in zip(specs, target):
        _register_collapse_alternates(fragment, chain, spec[1])
        name = _REGISTER_AA[letter]
        if not coot.mutate(fragment, *spec, name):
          raise RuntimeError(f"Mutation failed at {chain}:{spec[1]}.")
      if rigid_fit:
        fit_function = getattr(coot, "rigid_body_refine_by_atom_selection", None)
        if not callable(fit_function):
          raise RuntimeError("Rigid-body fitting is unavailable in this Coot build.")
        before_fit = _rigid_fit_atom_snapshot(fragment)
        fit_function(fragment, "//")
        coot.accept_regularizement()
        if _rigid_fit_atom_snapshot(fragment).keys() != before_fit.keys():
          raise RuntimeError("Rigid fitting changed fragment atom identities.")
      fit_rotamers()
      refinement = coot.refine_residues_py(fragment, specs)
      if refinement is False:
        raise RuntimeError("Coot could not refine this fragment.")
      coot.accept_regularizement()
      fit_rotamers()
      final_refinement = coot.refine_residues_py(fragment, specs)
      if final_refinement is False:
        raise RuntimeError("Coot could not perform the final fragment refinement.")
      coot.accept_regularizement()
    results = coot.map_to_model_correlation_per_residue_py(fragment, specs, 0, map_id) or []
    scores = {}
    for item in results:
      if not isinstance(item, (list, tuple)) or len(item) < 2:
        continue
      spec = _coot_residue_spec_from_spec(item[0])
      if spec is None:
        continue
      score = float(item[1])
      if math.isfinite(score) and list(spec) in specs:
        scores[tuple(spec)] = score
    if len(scores) != len(specs):
      raise RuntimeError(f"Only {len(scores)}/{len(specs)} residues could be scored; candidate excluded.")
    try:
      sidechain = _register_sidechain_correlation(fragment, specs, target, map_id)
    except Exception as error:
      print(f"Register {offset:+d}: sidechain CC unavailable: {error}")
      sidechain = {"mean": None, "count": 0, "eligible": sum(letter not in {"G", "A"} for letter in target)}
    coloured = False
    if _supports_direct_user_defined_colouring():
      rows = [(list(spec), _density_fit_score_to_colour_index(score)) for spec, score in scores.items()]
      coloured = bool(_apply_direct_user_defined_residue_colours(fragment, [], rows))
    return {"model": fragment, "offset": offset, "sequence": target,
            "mean": sum(scores.values())/len(scores), "minimum": min(scores.values()), "coloured": coloured,
            "sidechain": sidechain}
  except Exception:
    coot.close_molecule(fragment)
    raise
  finally:
    coot.set_refinement_immediate_replacement(previous_replacement)
    coot.set_imol_refinement_map(previous_map)


def compare_local_sequence_registers():
  active = _active_residue_or_status()
  if not active:
    return
  source, chain, centre, ins = active[:4]
  required = ("new_molecule_by_atom_selection", "renumber_residue_range", "mutate",
              "score_rotamers_py", "set_residue_to_rotamer_name", "refine_residues_py", "accept_regularizement",
              "map_to_model_correlation_per_residue_py", "set_imol_refinement_map",
              "refinement_immediate_replacement_state", "set_refinement_immediate_replacement",
              "set_mol_displayed", "mol_is_displayed", "close_molecule", "set_molecule_name",
              "delete_sequence_by_chain_id", "assign_sequence_from_string", "sequence_info_py")
  missing = [name for name in required if not callable(getattr(coot, name, None))]
  try:
    glib = _coot_gui_repository_module("GLib")
  except Exception:
    glib = None
  if missing or globals().get("Gtk") is None or glib is None:
    info_dialog("Register comparison requires the GTK interface and compatible editing APIs. " + ", ".join(missing))
    return
  if ins:
    info_dialog("Register comparison does not support an insertion-code centre.")
    return
  map_id = coot.imol_refinement_map()
  if map_id not in map_molecule_list() or map_is_difference_map(map_id):
    info_dialog("Choose a non-difference refinement map for register comparison.")
    return
  try:
    records = _register_chain_snapshot(source, chain)
    sequence, message = _register_reference_sequence(source, chain, records)
    origin = _register_sequence_origin(records, sequence) if sequence else None
    if message.startswith("Model-derived"):
      origin = min(r["number"] for r in records if not r["ins"] and r["name"] in _REGISTER_PROTEIN_NAMES - {"UNK"})
  except (ValueError, TypeError) as error:
    info_dialog(str(error))
    return
  centre_name = next((r["name"] for r in records if r["number"] == centre and r["ins"] == ins), "")
  needs_sequence_input = not sequence or origin is None or message.startswith("Model-derived")
  if not needs_sequence_input:
    try:
      _register_comparison_targets(records, centre, sequence, origin, 5, 9)
    except (ValueError, TypeError):
      needs_sequence_input = True
  centre_label = f"{chain}:{centre} {centre_name}".strip()
  window = Gtk.Window(title=f"Compare registers: #{source} | Centre {centre_label}")
  window.set_default_size(650, 480)
  box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6)
  for side in ("start", "end", "top", "bottom"):
    getattr(box, "set_margin_"+side)(8)
  label = Gtk.Label(label=f"Centre: {centre_label} | Model #{source}\n"
                         f"Region {chain}:{centre-4}-{centre+4} | Map #{map_id}")
  label.set_xalign(0)
  box.append(label)
  row = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
  origin_row = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
  origin_row.append(Gtk.Label(label="Sequence start"))
  origin_entry = Gtk.Entry()
  origin_entry.set_width_chars(7)
  origin_entry.set_max_width_chars(9)
  origin_entry.set_text(str(origin) if origin is not None else "")
  origin_row.append(origin_entry)
  row.append(Gtk.Label(label="Shift +/-"))
  shift_entry = Gtk.Entry()
  shift_entry.set_width_chars(3)
  shift_entry.set_max_width_chars(5)
  shift_entry.set_text("5")
  shift_entry.set_tooltip_text("Test every integer shift from -N to +N; 0 tests the original register only.")
  row.append(shift_entry)
  generate = Gtk.Button(label="Generate")
  row.append(generate)
  box.append(row)
  size_row = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
  size_row.append(Gtk.Label(label="Fragment size"))
  size_entry = Gtk.Entry()
  size_entry.set_width_chars(3)
  size_entry.set_max_width_chars(5)
  size_entry.set_text("9")
  size_entry.set_tooltip_text("Consecutive residues centred on the active residue. Even lengths extend one extra residue towards the C terminus.")
  size_row.append(size_entry)
  box.append(size_row)
  reference = Gtk.Expander(label="Reference sequence needs verification")
  reference.set_visible(needs_sequence_input)
  reference.set_expanded(needs_sequence_input)
  sequence_view = Gtk.TextView()
  sequence_view.set_wrap_mode(Gtk.WrapMode.WORD_CHAR)
  sequence_view.get_buffer().set_text(sequence)
  sequence_view.set_tooltip_text(message)
  sequence_scroll = Gtk.ScrolledWindow()
  sequence_scroll.set_min_content_height(70)
  sequence_scroll.set_child(sequence_view)
  sequence_controls = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6)
  sequence_controls.append(origin_row)
  sequence_controls.append(sequence_scroll)
  reference.set_child(sequence_controls)
  box.append(reference)
  rigid_fit_option = Gtk.CheckButton(label="Rigid fit before refinement")
  rigid_fit_available = callable(getattr(coot, "rigid_body_refine_by_atom_selection", None))
  rigid_fit_option.set_active(False)
  rigid_fit_option.set_sensitive(rigid_fit_available)
  rigid_fit_option.set_tooltip_text(
    "Mutate, rigid fit, select fixed-backbone rotamers, refine, select rotamers again, then final refine."
    if rigid_fit_available else "Rigid-body fitting is unavailable in this Coot build.")
  box.append(rigid_fit_option)
  progress = Gtk.Label(label="Isolated fragments; CC guides comparison, not certainty.")
  progress.set_wrap(True)
  progress.set_xalign(0)
  box.append(progress)
  headings = Gtk.Label(label="Shift   Sequence     Mean CC   Min CC   SC CC (n)")
  headings.set_xalign(0)
  headings.set_tooltip_text("SC CC: mean per-residue correlation beyond C-beta; coverage is scored/eligible residues. Gly/Ala excluded. Independent of contour.")
  headings.add_css_class("monospace")
  box.append(headings)
  rows = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=2)
  scroll = Gtk.ScrolledWindow()
  scroll.set_vexpand(True)
  scroll.set_child(rows)
  box.append(scroll)
  legend = Gtk.Label(label=_colour_band_legend_text("CC", 0.0, 1.0) if _supports_direct_user_defined_colouring()
                     else "Density-fit colouring unavailable in this build; scores still available.")
  legend.set_wrap(True)
  legend.set_xalign(0)
  legend_expander = Gtk.Expander(label="Density-fit colour scale")
  legend_expander.set_child(legend)
  box.append(legend_expander)
  buttons = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
  previous = Gtk.Button(label="Previous")
  next_button = Gtk.Button(label="Next")
  previous.set_sensitive(False)
  next_button.set_sensitive(False)
  buttons.append(previous)
  buttons.append(next_button)
  keep = Gtk.Button(label="Keep selected")
  keep.set_tooltip_text("Keep this fragment as a separate model; discard the other candidates.")
  keep.set_sensitive(False)
  discard = Gtk.Button(label="Discard all")
  buttons.append(keep)
  buttons.append(discard)
  box.append(buttons)
  state = {"closed": False, "candidates": [], "selected": None, "idle": None,
           "hidden_original": False, "index": 0, "targets": [], "snapshot": None, "buttons": [],
           "rigid_fit": False}
  def choose(candidate):
    if state["closed"] or not valid_model_molecule_qm(candidate["model"]):
      return
    for other in state["candidates"]:
      if valid_model_molecule_qm(other["model"]):
        coot.set_mol_displayed(other["model"], int(other is candidate))
    state["selected"] = candidate
    for button, other in state["buttons"]:
      if other is candidate:
        button.add_css_class("suggested-action")
      else:
        button.remove_css_class("suggested-action")
    keep.set_sensitive(state["idle"] is None)
    add_status_bar_text(f"Register {candidate['offset']:+d} | mean CC {candidate['mean']:.3f} | min CC {candidate['minimum']:.3f}")
  def navigate(direction):
    available = [c for c in state["candidates"] if valid_model_molecule_qm(c["model"])]
    if state["closed"] or not available:
      return
    index = next((i for i, c in enumerate(available) if c is state["selected"]), None)
    destination = (index + direction) % len(available) if index is not None else (0 if direction > 0 else -1)
    choose(available[destination])
  previous.connect("clicked", lambda *_: navigate(-1))
  next_button.connect("clicked", lambda *_: navigate(1))
  def finish(keep_selected=False):
    if state["closed"]:
      return
    state["closed"] = True
    if state["idle"] is not None:
      glib.source_remove(state["idle"])
      state["idle"] = None
    chosen = state["selected"] if keep_selected else None
    for candidate in state["candidates"]:
      if candidate is not chosen and valid_model_molecule_qm(candidate["model"]):
        _REGISTER_REFERENCE_ORIGINS.pop((candidate["model"], chain), None)
        coot.close_molecule(candidate["model"])
    if state["hidden_original"] and valid_model_molecule_qm(source):
      coot.set_mol_displayed(source, 1)
    window.destroy()
  def step():
    if state["closed"]:
      return False
    index = state["index"]
    offset, target = state["targets"][index]
    try:
      if map_id not in map_molecule_list() or _register_chain_snapshot(source, chain) != state["snapshot"]:
        raise ValueError("Source model or map changed; discard and restart the comparison.")
      candidate = _build_register_comparison_fragment(
        source, chain, state["start"], state["end"], offset, target, map_id,
        rigid_fit=state["rigid_fit"])
      try:
        _register_associate_reference(candidate["model"], chain, state["reference_sequence"], state["reference_origin"])
      except Exception:
        _REGISTER_REFERENCE_ORIGINS.pop((candidate["model"], chain), None)
        coot.close_molecule(candidate["model"])
        raise
      if state["closed"]:
        coot.close_molecule(candidate["model"])
        return False
      state["candidates"].append(candidate)
      button = Gtk.Button()
      sidechain = candidate["sidechain"]
      sidechain_text = f"{sidechain['mean']:.3f}" if sidechain["mean"] is not None else "n/a"
      candidate_label = Gtk.Label(label=f"{offset:+3d}     {target}    {candidate['mean']:.3f}    {candidate['minimum']:.3f}"
                                       f"    {sidechain_text} ({sidechain['count']}/{sidechain['eligible']})")
      candidate_label.set_xalign(0)
      candidate_label.add_css_class("monospace")
      button.set_child(candidate_label)
      button.set_tooltip_text(f"Centre {chain}:{centre+offset} | " +
                              ("Coloured by density fit" if candidate["coloured"] else "Density-fit colouring unavailable"))
      button.connect("clicked", lambda *_args, candidate=candidate: choose(candidate))
      rows.append(button)
      state["buttons"].append((button, candidate))
      previous.set_sensitive(True)
      next_button.set_sensitive(True)
      if state["selected"] is None:
        choose(candidate)
    except Exception as error:
      error_label = Gtk.Label(label=f"{offset:+d}  Failed: {error}")
      error_label.set_wrap(True)
      error_label.set_xalign(0)
      rows.append(error_label)
      print(f"Register comparison {offset:+d}: {error}")
    state["index"] += 1
    total = len(state["targets"])
    if state["index"] == total:
      state["idle"] = None
      keep.set_sensitive(state["selected"] is not None)
      available = [c for c in state["candidates"] if valid_model_molecule_qm(c["model"])]
      sidechain_available = [c for c in available if c["sidechain"]["mean"] is not None]
      best = max(sidechain_available or available,
                 key=lambda c: (c["sidechain"]["mean"], c["mean"]) if sidechain_available else (c["mean"],),
                 default=None)
      if best is not None:
        choose(best)
      metric_label = "SC CC" if sidechain_available else "mean CC (SC CC unavailable)"
      best_score = best["sidechain"]["mean"] if sidechain_available else (best["mean"] if best else None)
      progress.set_text(f"{len(state['candidates'])}/{total} built. Highest {metric_label}: {best['offset']:+d} ({best_score:.3f})." if best
                        else "No candidates succeeded; inspect the failed rows.")
      return False
    progress.set_text(f"Built {state['index']}/{total}; next shift {state['targets'][state['index']][0]:+d}")
    return True
  def begin(*_args, confirmed_size=None):
    if state["closed"] or state["idle"] is not None:
      return
    try:
      if map_id not in map_molecule_list():
        raise ValueError("The selected refinement map is no longer available.")
      pending_refinement = getattr(coot, "get_continue_updating_refinement_atoms_state", None)
      if callable(pending_refinement) and pending_refinement():
        raise ValueError("Accept or cancel the current refinement before comparing registers.")
      buffer = sequence_view.get_buffer()
      reference_sequence = buffer.get_text(buffer.get_start_iter(), buffer.get_end_iter(), False)
      snapshot = _register_chain_snapshot(source, chain)
      try:
        start, end, targets = _register_comparison_targets(
          snapshot, centre, reference_sequence, origin_entry.get_text(), shift_entry.get_text(), size_entry.get_text())
      except (ValueError, TypeError):
        reference.set_visible(True)
        reference.set_expanded(True)
        raise
    except (ValueError, TypeError) as error:
      progress.set_text(str(error))
      return
    fragment_size = end-start+1
    if fragment_size > 12 and confirmed_size != fragment_size:
      generate.set_sensitive(False)
      def cancel_generation():
        if not state["closed"]:
          generate.set_sensitive(True)
      def confirm_generation():
        if not state["closed"]:
          generate.set_sensitive(True)
          begin(confirmed_size=fragment_size)
      generic_confirm_dialog(
        "Generate large register fragments?",
        f"Fit {len(targets)} fragments of {fragment_size} residues each?\n"
        "Fragments longer than 12 residues may take a while to fit and refine.",
        "Cancel", cancel_generation, "Generate", confirm_generation,
        handle_close_function=cancel_generation,
      )
      return
    state.update(start=start, end=end, targets=targets, snapshot=snapshot,
                 reference_sequence=_register_clean_sequence(reference_sequence), reference_origin=int(origin_entry.get_text()),
                 rigid_fit=rigid_fit_available and bool(rigid_fit_option.get_active()))
    state["hidden_original"] = bool(coot.mol_is_displayed(source))
    if state["hidden_original"]:
      coot.set_mol_displayed(source, 0)
    generate.set_sensitive(False)
    origin_entry.set_sensitive(False)
    shift_entry.set_sensitive(False)
    size_entry.set_sensitive(False)
    label.set_text(f"Centre: {centre_label} | Model #{source}\nRegion {chain}:{start}-{end} | Map #{map_id}")
    sequence_view.set_editable(False)
    rigid_fit_option.set_sensitive(False)
    progress.set_text(f"Building register {targets[0][0]:+d}...")
    state["idle"] = glib.idle_add(step)
  generate.connect("clicked", begin)
  keep.connect("clicked", lambda *_: finish(True))
  discard.connect("clicked", lambda *_: finish(False))
  window.connect("close-request", lambda *_: finish(False) or True)
  window.set_child(box)
  window.present()


#Mutate range to poly-unk
def mutate_residue_range_by_click_a():
  def mutate_residue_range_by_click_b(res1,res2):
    clicked_range = _clicked_same_chain_residue_range(
      res1,
      res2,
      allow_single_residue=False,
      error_message="Start and end points must be in the same mol and chain!",
    )
    if clicked_range is None:
      return None
    mol_id, ch_id, res_start, res_end = clicked_range
    target_seq=(res_end-res_start+1)*"A"
    with _grouped_model_edit(mol_id):
      mutate_residue_range(mol_id,ch_id,res_start,res_end,target_seq)
      for resno in range(res_start,res_end+1):
        set_residue_name(mol_id,ch_id,resno,"","UNK")
  _script_user_defined_click(2, mutate_residue_range_by_click_b,
    prompt="Mutate range to UNK: click the start and end residues. Tab cancels.")

#Mutate range to polyala
def mutate_residue_range_by_click_ala_a():
  def mutate_residue_range_by_click_ala_b(res1,res2):
    clicked_range = _clicked_same_chain_residue_range(
      res1,
      res2,
      allow_single_residue=False,
      error_message="Start and end points must be in the same mol and chain!",
    )
    if clicked_range is None:
      return None
    mol_id, ch_id, res_start, res_end = clicked_range
    target_seq=(res_end-res_start+1)*"A"
    with _grouped_model_edit(mol_id):
      mutate_residue_range(mol_id,ch_id,res_start,res_end,target_seq)
      for resno in range(res_start,res_end+1):
        set_residue_name(mol_id,ch_id,resno,"","ALA")
  _script_user_defined_click(2, mutate_residue_range_by_click_ala_b,
    prompt="Mutate range to ALA: click the start and end residues. Tab cancels.")
  
#Force addition of residue - useful when
#Coot says "No acceptable position found"
# but density is clear.
def force_add_terminal_residue():
  def force_addition(res1):
    mol_id=_click_spec_imol(res1)
    ch_id=_click_spec_chain_id(res1)
    res_no=_click_spec_res_no(res1)
    ins_code=_click_spec_ins_code(res1)
    if mol_id < 0 or ins_code:
      info_dialog("Select a uniquely identified terminal residue without an insertion code.")
      return 0
    return force_add_terminal_residue_noclick_phi_psi(mol_id, ch_id, res_no, -57.82, -47)
  _script_user_defined_click(1,force_addition)
  
def force_add_terminal_residue_noclick(mol_id,ch_id,res_no):
  return force_add_terminal_residue_noclick_phi_psi(mol_id, ch_id, res_no, -57.82, -47)

def force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,res_no,phi,psi,res_type="auto"):
  candidates = (res_no-1, res_no+1)
  before = {number for number in candidates if residue_exists_qm(mol_id, ch_id, number, "")}
  add_terminal_residue_using_phi_psi(mol_id,ch_id,res_no,res_type,float(phi),float(psi))
  created = [number for number in candidates if number not in before
             and residue_exists_qm(mol_id, ch_id, number, "")]
  for number in created:
    set_b_factor_residue_range(mol_id, ch_id, number, number, default_new_atoms_b_factor())
  sort_residues(mol_id)
  return int(len(created) == 1)

def force_add_terminal_residue_noclick_strand(mol_id,ch_id,res_no):
  return force_add_terminal_residue_noclick_phi_psi(mol_id, ch_id, res_no, -139, 135)

#Paul
def key_binding_terminal_spin():
    active_atom = active_residue()
    if not active_atom:
       add_status_bar_text("No active atom")
       return None
    else:
       imol = active_atom[0]
       residue_spec = _coot_residue_spec_from_spec(atom_spec_to_residue_spec(active_atom))
       if residue_spec is None:
         add_status_bar_text("Could not resolve the active residue for terminal spin")
         return None
       print(('spin_N {} {} {}'.format(imol, residue_spec, 120)))
       spin_N_py(imol, list(residue_spec), 120)

#This works, but should alter to be more flexible. ideally, have one key to tweak phi, one to tweak psi.
#Will need two global cycle variables - residue_phi_cycle and residue_psi_cycle - as well as two variables to keep the current value of phi and psi.
residue_phi_cycle=0
residue_psi_cycle=0
current_phi=-60
current_psi=-50
RESIDUE_TORSION_CYCLE_VALUES = (-180, -150, -120, -90, -60, -30, 0, 30, 60, 90, 120)
#should chamnge this to incorporate measurement of starting phi/psi using get_torsion (and maybe setting using set_torsion?)
#get_torsion(0,["A",2393,""," C  ",""], ["A",2394,"", " N  ", ""], ["A", 2394, "", " CA ", ""], ["A", 2394, "", " C  ",""])
#set_torsion(imol, chain_id, res_no, ins_code_alt_conf, atom_name_1,atom_name_2, atom_name_3, atom_name_4)
def _cycle_terminal_residue_torsion(vary_phi):
  global residue_phi_cycle
  global residue_psi_cycle
  global current_phi
  global current_psi
  growth_context = _active_terminal_growth_context_or_status()
  if not growth_context:
    return None
  cycle_index = residue_phi_cycle if vary_phi else residue_psi_cycle
  chosen_angle = RESIDUE_TORSION_CYCLE_VALUES[cycle_index]
  mol_id = growth_context["mol_id"]
  ch_id = growth_context["chain_id"]
  atom_name = growth_context["atom_name"]
  anchor_resno = growth_context["first_in_seg"] if growth_context["grow_from_n_term"] else growth_context["last_in_seg"]
  insertion_resno = anchor_resno + 1 if growth_context["grow_from_n_term"] else anchor_resno - 1
  specs = growth_context["segment_specs"]
  if len(specs) < 2 or any(spec[2] for spec in specs):
    info_dialog("Terminal torsion cycling needs a neighbouring anchor and no insertion codes.")
    return 0
  endpoint = (ch_id, anchor_resno, "")
  neighbour = specs[1] if growth_context["grow_from_n_term"] else specs[-2]
  if neighbour[1] != insertion_resno:
    info_dialog("Terminal torsion cycling needs consecutive endpoint numbering.")
    return 0
  original_atoms = residue_info_py(mol_id, *endpoint) or []
  if any(atom[0][1] for atom in original_atoms):
    info_dialog("Terminal torsion cycling cannot safely rebuild alternate conformers.")
    return 0
  original_type = residue_name(mol_id, *endpoint)
  next_cycle_index = (cycle_index + 1) % len(RESIDUE_TORSION_CYCLE_VALUES)
  if vary_phi:
    phi = chosen_angle
    psi = current_psi
    residue_phi_cycle = next_cycle_index
  else:
    phi = current_phi
    psi = chosen_angle
    residue_psi_cycle = next_cycle_index
  fragment = None
  try:
    fragment = _copy_residue_specs_verified(mol_id, specs)
    set_mol_displayed(fragment, 0)
    delete_residue(fragment, *endpoint)
    if not force_add_terminal_residue_noclick_phi_psi(
        fragment, ch_id, insertion_resno, phi, psi, original_type):
      raise RuntimeError("Coot could not rebuild the terminal residue.")
    if not residue_exists_qm(fragment, *endpoint) or residue_name(fragment, *endpoint) != original_type:
      raise RuntimeError("The rebuilt endpoint has an unexpected identity.")
    with _grouped_model_edit(mol_id):
      if not coot.replace_fragment(mol_id, fragment, _residue_spec_to_cid_compat(endpoint)):
        raise RuntimeError("Coot rejected the rebuilt endpoint; use Undo if necessary.")
  except Exception as error:
    info_dialog(f"Terminal torsion cycling stopped: {error}")
    return 0
  finally:
    if fragment is not None and valid_model_molecule_qm(fragment):
      close_molecule(fragment)
  set_go_to_atom_molecule(mol_id)
  set_go_to_atom_chain_residue_atom_name(ch_id, anchor_resno, atom_name)
  current_phi = phi
  current_psi = psi
  return 1


def cycle_residue_phi():
  return _cycle_terminal_residue_torsion(vary_phi=True)


def cycle_residue_psi():
  return _cycle_terminal_residue_torsion(vary_phi=False)


def _mutate_all_residue_type(source_resname, target_resname):
  mol_id = _active_molecule_or_status()
  if mol_id is None:
    return None
  map_id = imol_refinement_map()
  with _grouped_model_edit(mol_id):
    for ch_id in chain_ids(mol_id):
      for resn in range(0, chain_n_residues(ch_id, mol_id)):
        if resname_from_serial_number(mol_id, ch_id, resn) == source_resname:
          seqnum=seqnum_from_serial_number(mol_id,ch_id,resn)
          ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,resn))
          delete_residue_with_full_spec(mol_id,1,ch_id,seqnum,ins_id,"B")
          delete_residue_with_full_spec(mol_id,1,ch_id,seqnum,ins_id,"C")
          mutate(mol_id,ch_id,seqnum,ins_id,target_resname)
          if (map_id!=-1):
            coot.auto_fit_best_rotamer(mol_id,ch_id,seqnum,ins_id,"",map_id,1,0.01)


#Mutate mets to MSE
def mutate_all_mets_to_mse():
  return _mutate_all_residue_type("MET", "MSE")
  
#Mutate MSEs to MET
def mutate_all_mse_to_met():
  return _mutate_all_residue_type("MSE", "MET")
  
#Shorten loop by one residue
def shorten_loop():
  active_atom=_active_residue_or_status()
  if not active_atom:
    return None
  mol_id=active_atom[0]
  ch_id=active_atom[1]
  resn=active_atom[2]
  ins_code=active_atom[3]
  delete_residue(mol_id,ch_id,resn,ins_code)
  first_res=seqnum_from_serial_number(mol_id,"%s"%(ch_id),0)
  if not _renumber_residue_range_or_dialog(
    mol_id,
    ch_id,
    first_res,
    resn,
    1,
    "Failed to shorten the loop because the renumbering step was not possible.\n\n"
    "This usually means the requested change would create overlapping sequence numbering.",
  ):
    return None
  _refresh_extra_restraints_display(mol_id)
  r1=resn-1
  r2=resn+2
  previous_replacement = refinement_immediate_replacement_state()
  try:
    set_refinement_immediate_replacement(1)
    refine_zone(mol_id,ch_id,r1,r2,"")
    accept_regularizement()
  finally:
    set_refinement_immediate_replacement(previous_replacement)

#Lengthen loop by one residue
def lengthen_loop():
  active_atom=_active_residue_or_status()
  if not active_atom:
    return None
  mol_id=active_atom[0]
  ch_id=active_atom[1]
  resn=active_atom[2]
  first_res=seqnum_from_serial_number(mol_id,"%s"%(ch_id),0)
  if not _renumber_residue_range_or_dialog(
    mol_id,
    ch_id,
    first_res,
    resn,
    -1,
    "Failed to lengthen the loop because the renumbering step was not possible.\n\n"
    "This usually means the requested change would create overlapping sequence numbering.",
  ):
    return None
  sort_residues(mol_id)
  _refresh_extra_restraints_display(mol_id)
  r1=resn-1
  r2=resn
  previous_replacement = refinement_immediate_replacement_state()
  set_refinement_immediate_replacement(1)
  try:
    status = add_residue_by_map_fit(mol_id,ch_id,r1,"auto",1)
    if not status:
      info_dialog("Failed to add a residue into the loop gap")
      return None
    refine_zone(mol_id,ch_id,r1,r2+1,"")
    accept_regularizement()
  finally:
    set_refinement_immediate_replacement(previous_replacement)

#Get fractional coordinates of active atom. Useful when inspecting heavy atom sites.
def get_fract_coords():
  a=_active_residue_or_status()
  if not a:
    return None
  x_cart=atom_specs(a[0],a[1],a[2],a[3],a[4],a[5])[3]
  y_cart=atom_specs(a[0],a[1],a[2],a[3],a[4],a[5])[4]
  z_cart=atom_specs(a[0],a[1],a[2],a[3],a[4],a[5])[5]
  mol_id=a[0]
  molecule_cell = cell(mol_id)
  cell_a=molecule_cell[0]
  cell_b=molecule_cell[1]
  cell_c=molecule_cell[2]
  cell_alpha=math.radians(molecule_cell[3])
  cell_beta=math.radians(molecule_cell[4])
  cell_gamma=math.radians(molecule_cell[5])
  cos_alpha_star=(math.cos(cell_beta)*math.cos(cell_gamma)-math.cos(cell_alpha))/(math.sin(cell_beta)*math.sin(cell_gamma))
  sin_alpha_star=math.sqrt(1-cos_alpha_star**2)
  z_fract=z_cart/(cell_c*math.sin(cell_beta)*sin_alpha_star)
  y_fract=(y_cart-(-1*cell_c*math.sin(cell_beta)*cos_alpha_star)*z_fract)/(cell_b*math.sin(cell_gamma))
  x_fract=(x_cart-(cell_b*math.cos(cell_gamma))*y_fract-(cell_c*math.cos(cell_beta))*z_fract)/cell_a
  x_y_z_string=str("("+str("%.3f" % x_fract)+","+str("%.3f" % y_fract)+","+str("%.3f" % z_fract)+")")
  info_dialog("The fractional coordinates of the active atom are: %s"%(x_y_z_string))
  
#Go to center of scroll wheel map.
def goto_center_of_map():
  map_id = _scrollable_map_or_status()
  if map_id is None:
    info_dialog("You need a map!")
    return None
  a=float(cell(map_id)[0])
  b=float(cell(map_id)[1])
  c=float(cell(map_id)[2])
  x=a*0.5
  y=b*0.5
  z=c*0.5
  set_rotation_centre(x,y,z)
    
#Switch all models to CA representation
def all_mols_to_ca():
  for mol_id in molecule_number_list():
    graphics_to_ca_plus_ligands_representation(mol_id)
    
#Set b-factor color scaling based on mean B of active mol
def autoscale_b_factor():
  mol_id = _active_molecule_or_status()
  if mol_id is None:
    return None
  mean_b=average_temperature_factor(mol_id)
  scale_fac=50/mean_b
  set_b_factor_bonds_scale_factor(mol_id,scale_fac)
  graphics_to_b_factor_representation(mol_id)
    
#Color molecule by rotamer and missing atom outliers
def _iter_rotamer_outlier_colours(mol_id, defer_apply):
  missing_atoms_list=[]
  missing_atoms_colour=2
  rotamer_outlier_list=[]
  rotamer_outlier_colour=34
  blank_colour=0
  for x in missing_atom_info(mol_id):
    missing_atoms_list.append((x, missing_atoms_colour))
  specs = _all_residue_specs_for_colouring(mol_id)
  for index, (ch_id, resn, ins) in enumerate(specs):
    rot_prob = rotamer_score(mol_id, ch_id, resn, ins, "")
    colour = rotamer_outlier_colour if 0.0 < rot_prob < 0.5 else blank_colour
    rotamer_outlier_list.append(([ch_id, resn, ins], colour))
    if (index + 1) % 50 == 0:
      yield f"Scoring rotamers: {index+1}/{len(specs)}", None, True
  colour_list = rotamer_outlier_list + missing_atoms_list
  result = _apply_user_defined_residue_colours(mol_id, [], colour_list, defer_apply=defer_apply)
  if defer_apply:
    return (yield from result)
  return result


def color_rotamer_outliers_and_missing_atoms(mol_id, defer_apply=False):
  steps = _iter_rotamer_outlier_colours(mol_id, defer_apply)
  return steps if defer_apply else _consume_colour_steps(steps)


def _residue_names_for_colouring(mol_id):
  try:
    return [(list(spec), record["name"]) for spec, record in _model_coordinate_snapshot(mol_id).items()]
  except ValueError:
    return [(spec, residue_name(mol_id, *spec)) for spec in _all_residue_specs_for_colouring(mol_id)]


def color_polars_and_hphobs(mol_id, defer_apply=False):
  hphob_list=["CYS","ILE","LEU","VAL","TYR","MET","PHE","TRP","ALA"]
  polar_list=["SER","ASN","GLN","HIS","ARG","LYS","GLU","ASP","THR"]
  #based these on Moon&Fleming PNAS 2011 and MacCallum TIBS 2011
  #Gly/Pro colored differently because they are conformationally "special" residues
  polar_colour=5 #light blue
  hphob_colour=28 #orange
  gly_colour=34 #magenta
  pro_colour=15 #green
  blank_colour=0 #light gray
  polar_res_list=[]
  hphob_res_list=[]
  gly_res_list=[]
  pro_res_list=[]
  blank_res_list=[]
  for spec, resname_here in _residue_names_for_colouring(mol_id):
    if resname_here in polar_list:
      polar_res_list.append((spec, polar_colour))
    elif resname_here in hphob_list:
      hphob_res_list.append((spec, hphob_colour))
    elif resname_here == "GLY":
      gly_res_list.append((spec, gly_colour))
    elif resname_here == "PRO":
      pro_res_list.append((spec, pro_colour))
    else:
      blank_res_list.append((spec, blank_colour))
  colour_list = polar_res_list + hphob_res_list + gly_res_list + pro_res_list
  return _apply_user_defined_residue_colours(mol_id, blank_res_list, colour_list, defer_apply=defer_apply)
  
def color_by_charge(mol_id, defer_apply=False):
  pos_list=["ARG","LYS","HIS"]
  neg_list=["GLU","ASP"]
  pos_colour=4
  neg_colour=31
  blank_colour=0
  pos_res_list=[]
  neg_res_list=[]
  blank_res_list=[]
  for spec, resname_here in _residue_names_for_colouring(mol_id):
    if resname_here in pos_list:
      pos_res_list.append((spec, pos_colour))
    elif resname_here in neg_list:
      neg_res_list.append((spec, neg_colour))
    else:
      blank_res_list.append((spec, blank_colour))
  colour_list = pos_res_list + neg_res_list
  return _apply_user_defined_residue_colours(mol_id, blank_res_list, colour_list, defer_apply=defer_apply)


def _append_chain_residue_colours(colour_rows, mol_id, ch_id, colour_index):
  if _supports_direct_user_defined_colouring():
    colour_rows.append((f"//{ch_id}", colour_index))
    return
  for sn in _iter_chain_serial_numbers(mol_id, ch_id):
    resn = seqnum_from_serial_number(mol_id, ch_id, sn)
    ins_id = str(insertion_code_from_serial_number(mol_id, ch_id, sn))
    colour_rows.append(([ch_id, resn, ins_id], colour_index))


def uncolor_other_chains():
  chain_context = _active_chain_context_or_status()
  if not chain_context:
    return None
  mol_id = chain_context["mol_id"]
  ch_id_here = chain_context["chain_id"]
  blank_colour=0
  blank_res_list=[]
  for ch_id in chain_ids(mol_id):
    if ch_id!=ch_id_here:
      _append_chain_residue_colours(blank_res_list, mol_id, ch_id, blank_colour)
  if _supports_direct_user_defined_colouring():
    _set_user_defined_atom_colour_by_residue_rows_py(mol_id, blank_res_list)
    _graphics_to_legacy_user_defined_representation(mol_id, blank_res_list)
  else:
    _apply_user_defined_residue_colours(mol_id, blank_res_list, [])

def color_active_chain():
  chain_context = _active_chain_context_or_status()
  if not chain_context:
    return None
  mol_id = chain_context["mol_id"]
  ch_id_here = chain_context["chain_id"]
  blank_colour=0
  chain_colour=22 #yellow
  blank_res_list=[]
  chain_res_list=[]
  for ch_id in chain_ids(mol_id):
    if ch_id!=ch_id_here:
      _append_chain_residue_colours(blank_res_list, mol_id, ch_id, blank_colour)
    else:
      _append_chain_residue_colours(chain_res_list, mol_id, ch_id, chain_colour)
  _apply_user_defined_residue_colours(mol_id, blank_res_list, chain_res_list)

def color_active_chain_by_num(chain_colour):
  chain_context = _active_chain_context_or_status()
  if not chain_context:
    return None
  mol_id = chain_context["mol_id"]
  ch_id_here = chain_context["chain_id"]
  blank_colour=0
  blank_res_list=[]
  chain_res_list=[]
  for ch_id in chain_ids(mol_id):
    if ch_id!=ch_id_here:
      _append_chain_residue_colours(blank_res_list, mol_id, ch_id, blank_colour)
    else:
      _append_chain_residue_colours(chain_res_list, mol_id, ch_id, chain_colour)
  _apply_user_defined_residue_colours(mol_id, blank_res_list, chain_res_list)
    
def color_protein_na(mol_id, defer_apply=False):
  blank_colour=0
  protein_colour=22 #yellow
  na_colour=31
  protein_res_list=[]
  na_res_list=[]
  blank_res_list=[]
  for ch_id in chain_ids(mol_id):
    if is_nucleotide_chain_p(mol_id,ch_id):
      _append_chain_residue_colours(na_res_list, mol_id, ch_id, na_colour)
    elif is_protein_chain_p(mol_id,ch_id):
      _append_chain_residue_colours(protein_res_list, mol_id, ch_id, protein_colour)
    else:
      _append_chain_residue_colours(blank_res_list, mol_id, ch_id, blank_colour)
  colour_list = protein_res_list + na_res_list
  return _apply_user_defined_residue_colours(mol_id, blank_res_list, colour_list, defer_apply=defer_apply)


def color_waters(mol_id, defer_apply=False):
  water_colour=31
  blank_colour=0
  water_list=[]
  blank_res_list=[]
  for spec, resname_here in _residue_names_for_colouring(mol_id):
    colour_entry = (spec, water_colour if resname_here=="HOH" else blank_colour)
    if resname_here=="HOH":
      water_list.append(colour_entry)
    else:
      blank_res_list.append(colour_entry)
  return _apply_user_defined_residue_colours(mol_id, blank_res_list, water_list, defer_apply=defer_apply)


def color_rotamer_outliers_and_missing_atoms_for_active_molecule():
  mol_id = _active_molecule_or_status()
  if mol_id is None:
    return None
  return _start_user_defined_colour_action(mol_id, color_rotamer_outliers_and_missing_atoms, "Rotamer / missing-atom colouring")


def color_polars_and_hphobs_for_active_molecule():
  mol_id = _active_molecule_or_status()
  if mol_id is None:
    return None
  return _start_user_defined_colour_action(mol_id, color_polars_and_hphobs, "Polarity colouring")


def color_by_charge_for_active_molecule():
  mol_id = _active_molecule_or_status()
  if mol_id is None:
    return None
  return _start_user_defined_colour_action(mol_id, color_by_charge, "Charge colouring")


def color_protein_na_for_active_molecule():
  mol_id = _active_molecule_or_status()
  if mol_id is None:
    return None
  return _start_user_defined_colour_action(mol_id, color_protein_na, "Protein / nucleic-acid colouring")


def color_waters_for_active_molecule():
  mol_id = _active_molecule_or_status()
  if mol_id is None:
    return None
  return _start_user_defined_colour_action(mol_id, color_waters, "Water colouring")


CHAIN_BREAK_PROTEIN_RESNAMES = {'ALA','UNK','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','ILE','LEU','LYS','MET','MSE','PHE','PRO','SER','THR','TRP','TYR','VAL'}
CHAIN_BREAK_NA_RESNAMES = {'A','C','T','G','U'}


def _chain_break_config_for_resname(resname):
  if resname in CHAIN_BREAK_PROTEIN_RESNAMES:
    return ("CA", 3.8)
  if resname in CHAIN_BREAK_NA_RESNAMES:
    return ("P", 5.9)
  return None


def _atom_xyz_or_none(mol_id, ch_id, resno, atom_name):
  try:
    atom_spec = atom_specs(mol_id, ch_id, resno, "", atom_name, "")
    return (atom_spec[-3], atom_spec[-2], atom_spec[-1])
  except Exception:
    return None


def _chain_break_gap_color(res_missing, distance_per_residue, warning_threshold):
  if res_missing <= 15:
    gap_color = "gray"
  elif res_missing < 50:
    gap_color = "orange"
  else:
    gap_color = "red"
  if distance_per_residue > warning_threshold:
    gap_color = "cyan"
  return gap_color


def _chain_break_label(ch_id, start_resno, end_resno, res_missing, distance, distance_per_residue, warning_threshold):
  label_string = "{ch_id}  {res_start}...{res_end} ({res_missing} missing residues)".format(
    ch_id=ch_id,
    res_start=start_resno,
    res_end=end_resno,
    res_missing=res_missing,
  )
  if distance_per_residue > warning_threshold:
    label_string = "MIND THE GAP! {ch_id}  {res_start}...{res_end} ({res_missing} missing residues for {distance:6.1f} A)".format(
      ch_id=ch_id,
      res_start=start_resno,
      res_end=end_resno,
      res_missing=res_missing,
      distance=distance,
    )
  if res_missing >= 50:
    label_string = "!!!  {ch_id}  {res_start}...{res_end} ({res_missing} missing residues)".format(
      ch_id=ch_id,
      res_start=start_resno,
      res_end=end_resno,
      res_missing=res_missing,
    )
  return label_string


def _highlight_chain_breaks_for_molecule(mol_id, collect_missing_segments):
  clear_ball_and_stick(mol_id)
  missing_segments_list = []
  snapshot = _model_coordinate_snapshot(mol_id)
  with _grouped_model_edit(mol_id):
    obj_number = generic_object_with_name("chain_breaks_{mol_id}".format(mol_id=mol_id))
    generic_object_clear(obj_number)
    for ch_id in chain_ids(mol_id):
      rows = _segment_chain_snapshot(mol_id, ch_id, snapshot)
      breaks = [left["name"] in _SEGMENT_POLYMER_NAMES and right["name"] in _SEGMENT_POLYMER_NAMES
                and not _segment_rows_connected(left, right) for left, right in zip(rows, rows[1:])]
      previous_terminus_xyz = {}
      for resn, row in enumerate(rows):
        resname = row["name"]
        config = _chain_break_config_for_resname(resname)
        if not config:
          continue
        atom_name, warning_threshold = config
        seqnum = row["number"]
        if seqnum == -10000:
          continue
        if resn == 0 or resn == len(rows)-1 or rows[resn-1]["name"] not in _SEGMENT_POLYMER_NAMES or (resn+1 < len(rows) and rows[resn+1]["name"] not in _SEGMENT_POLYMER_NAMES):
          sel_string = _residue_spec_to_cid_compat((ch_id, seqnum, row["ins"])) + "/" + atom_name
          make_ball_and_stick(mol_id, sel_string, 0, 0.5, 1)
        if resn < len(breaks) and breaks[resn]:
          previous_xyz = row["xyz"].get(atom_name.strip())
          if previous_xyz is not None:
            previous_terminus_xyz[atom_name] = previous_xyz
        if resn and breaks[resn-1]:
          previous_xyz = previous_terminus_xyz.get(atom_name)
          current_xyz = row["xyz"].get(atom_name.strip())
          previous_seqnum = rows[resn-1]["number"]
          if previous_xyz is None or current_xyz is None or previous_seqnum == -10000:
            continue
          res_missing = seqnum - previous_seqnum
          if res_missing <= 0:
            continue
          x_here, y_here, z_here = previous_xyz
          x_mn, y_mn, z_mn = current_xyz
          x_mid=(x_mn+x_here)/2
          y_mid=(y_mn+y_here)/2
          z_mid=(z_mn+z_here)/2
          distance=((x_here-x_mn)**2+(y_here-y_mn)**2+(z_here-z_mn)**2)**0.5
          distance_per_residue=distance/res_missing
          gap_color = _chain_break_gap_color(res_missing, distance_per_residue, warning_threshold)
          if collect_missing_segments:
            label_string = _chain_break_label(
              ch_id,
              previous_seqnum,
              seqnum,
              res_missing,
              distance,
              distance_per_residue,
              warning_threshold,
            )
            missing_segments_list.append([label_string, x_mid, y_mid, z_mid])
          to_generic_object_add_dashed_line(obj_number, gap_color, 4, 3, x_here, y_here, z_here, x_mn, y_mn, z_mn)
    set_display_generic_object(obj_number,1)
    try:
      attach_generic_object_to_molecule(obj_number, mol_id)
    except NameError:
      info_dialog("attach_generic_object_to_molecule is only present in Coot r6057 and later, sorry.")
  return missing_segments_list


#Highlight various items with ball-and-stick and lines
def highlight_chain_breaks():
  mol_id = _active_molecule_or_status()
  if mol_id is None:
    return None
  missing_segments_list = _highlight_chain_breaks_for_molecule(mol_id, True)
  interesting_things_gui("Missing segments",missing_segments_list)


def highlight_all_chain_breaks():
  for mol_id in model_molecule_list():
    _highlight_chain_breaks_for_molecule(mol_id, False)

#Place helix and add helix restraints
def place_helix_with_restraints():
  place_helix_here()

#Merge two chains (throw error msg if they overlap in numbering)
#This needs fixing - says selections overlap in cases where sel1 has one segment on either side from sel2
#Fixed now, I think.
def merge_chains():
  def merge_chains_by_click(res1,res2):
    mol_id_1=_click_spec_imol(res1)
    mol_id_2=_click_spec_imol(res2)
    ch_id_1=_click_spec_chain_id(res1)
    ch_id_2=_click_spec_chain_id(res2)
    if mol_id_1 < 0 or mol_id_1 != mol_id_2 or ch_id_1 == ch_id_2:
      info_dialog("Select two different chains in the same model.")
      return
    first_res_sel1=first_residue(mol_id_1,ch_id_1)
    last_res_sel1=last_residue(mol_id_1,ch_id_1)
    first_res_sel2=first_residue(mol_id_2,ch_id_2)
    last_res_sel2=last_residue(mol_id_2,ch_id_2)
    size_2=abs(last_res_sel2-first_res_sel2)
    size_1=abs(last_res_sel1-first_res_sel1)
    if (mol_id_1!=mol_id_2) or (ch_id_1==ch_id_2):
      info_dialog("No can do, chains must be in the same mol and have non-overlapping ranges!")
    elif (mol_id_1==mol_id_2) and (ch_id_1!=ch_id_2) and (size_1<=size_2):
      out=change_chain_id_with_result(mol_id_1,ch_id_1,ch_id_2,1,first_res_sel1,last_res_sel1)
      if out[0]==0:
        info_dialog(out[1])
    elif (mol_id_1==mol_id_2) and (ch_id_1!=ch_id_2) and (size_2<size_1):
      out=change_chain_id_with_result(mol_id_1,ch_id_2,ch_id_1,1,first_res_sel2,last_res_sel2)
      if out[0]==0:
        info_dialog(out[1])
    else:
      info_dialog("No can do, chains must be in the same mol and have non-overlapping ranges!")
  _script_user_defined_click(2, merge_chains_by_click,
    prompt="Merge chains: click a residue in each chain. Tab cancels.")

#Fit all segments to map
def rigid_body_fit_segments():
  if (imol_refinement_map()==-1):
    info_dialog("You must set a refinement map!")
    return None
  mol_id = _active_molecule_or_status()
  if mol_id is None:
    return None
  previous_immediate_replacement = refinement_immediate_replacement_state()
  with _grouped_model_edit(mol_id):
    try:
      for specs in _model_segment_specs(mol_id):
        set_refinement_immediate_replacement(1)
        _rigid_fit_residue_specs(mol_id, specs)
        accept_regularizement()
    finally:
      set_refinement_immediate_replacement(previous_immediate_replacement)
      
#Fit current segment
def fit_this_segment():
  if (imol_refinement_map()==-1):
    info_dialog("You must set a refinement map!")
    return None
  segment_context = _active_segment_context_or_status()
  if not segment_context:
    return None
  mol_id = segment_context["mol_id"]
  previous_immediate_replacement = refinement_immediate_replacement_state()
  with _grouped_model_edit(mol_id):
    try:
      set_refinement_immediate_replacement(1)
      _rigid_fit_residue_specs(mol_id, segment_context["segment_specs"])
      accept_regularizement()
    finally:
      set_refinement_immediate_replacement(previous_immediate_replacement)


def _rigid_fit_residue_specs(mol_id, specs):
  fit = getattr(coot, "rigid_body_refine_by_atom_selection", None)
  if callable(fit):
    return fit(mol_id, _selection_for_residue_specs(specs))
  chain = specs[0][0]
  start, end = min(spec[1] for spec in specs), max(spec[1] for spec in specs)
  actual = {tuple(spec) for spec in _tq_chain_residue_specs(mol_id, chain) if start <= spec[1] <= end}
  if any(spec[2] for spec in specs) or actual != set(specs):
    raise ValueError("This Coot build cannot rigid-fit this insertion-coded segment safely.")
  return rigid_body_refine_zone(mol_id, chain, start, end)

#Set default b-fac for new atoms to mean B for active mol
def set_new_atom_b_fac_to_mean():
  mol_id = _active_molecule_or_status()
  if mol_id is None:
    return None
  mean_b=average_temperature_factor(mol_id)
  set_default_temperature_factor_for_new_atoms(mean_b)
  

def _active_terminal_growth_context_or_status():
  active_atom = _active_atom_context_or_status()
  if not active_atom:
    return None
  mol_id = active_atom["mol_id"]
  ch_id = active_atom["chain_id"]
  resn = active_atom["resno"]
  specs = _segment_specs_for_residue(mol_id, ch_id, resn, active_atom["ins_code"])
  if not specs:
    add_status_bar_text("No connected polymer segment at the active residue")
    return None
  first_in_seg = specs[0][1]
  last_in_seg = specs[-1][1]
  active_atom["segment_specs"] = specs
  active_atom["first_in_seg"] = first_in_seg
  active_atom["last_in_seg"] = last_in_seg
  active_atom["grow_from_n_term"] = abs(first_in_seg - resn) <= abs(last_in_seg - resn)
  return active_atom


def _run_terminal_growth_action(builder):
  growth_context = _active_terminal_growth_context_or_status()
  if not growth_context:
    return None
  anchor = (growth_context["segment_specs"][0] if growth_context["grow_from_n_term"]
            else growth_context["segment_specs"][-1])
  if anchor[2]:
    info_dialog("Terminal growth at an insertion-coded residue is not supported by Coot's numbered builder.")
    return None
  set_go_to_atom_molecule(growth_context["mol_id"])
  set_new_atom_b_fac_to_mean()
  if growth_context["grow_from_n_term"]:
    anchor_resno = growth_context["first_in_seg"]
    destination_resno = anchor_resno - 1
  else:
    anchor_resno = growth_context["last_in_seg"]
    destination_resno = anchor_resno + 1
  set_go_to_atom_chain_residue_atom_name(growth_context["chain_id"], anchor_resno, growth_context["atom_name"])
  builder(growth_context, anchor_resno)
  sort_residues(growth_context["mol_id"])
  set_go_to_atom_chain_residue_atom_name(growth_context["chain_id"], destination_resno, growth_context["atom_name"])
  return 1
  
  
#Shortcut for adding terminal residue (to bind to key)
def add_term_shortcut():
  if imol_refinement_map()==-1:
    info_dialog("You must set a refinement map!")
    return None
  return _run_terminal_growth_action(
    lambda growth_context, anchor_resno: add_terminal_residue(
      growth_context["mol_id"],
      growth_context["chain_id"],
      anchor_resno,
      _register_terminal_residue_type(
        growth_context["mol_id"], growth_context["chain_id"],
        anchor_resno + (-1 if growth_context["grow_from_n_term"] else 1)),
      1,
    )
  )
    
def add_term_shortcut_force():
  return _run_terminal_growth_action(
    lambda growth_context, anchor_resno: force_add_terminal_residue_noclick(
      growth_context["mol_id"],
      growth_context["chain_id"],
      anchor_resno,
    )
  )
    
def add_term_shortcut_force_strand():
  return _run_terminal_growth_action(
    lambda growth_context, anchor_resno: force_add_terminal_residue_noclick_strand(
      growth_context["mol_id"],
      growth_context["chain_id"],
      anchor_resno,
    )
  )

#Fix all atoms in active residue
def fix_active_residue():
  active_atom = _active_atom_context_or_status()
  if not active_atom:
    return None
  residue_info_list=residue_info_py(
    active_atom["mol_id"],
    active_atom["chain_id"],
    active_atom["resno"],
    active_atom["ins_code"],
  )
  atom_list=[]
  for item in residue_info_list:
    atom_name=item[0][0]
    atom_info=[
      active_atom["chain_id"],
      active_atom["resno"],
      active_atom["ins_code"],
      atom_name,
      active_atom["alt_conf"],
    ]
    atom_list.append(atom_info)
  mark_multiple_atoms_as_fixed(active_atom["mol_id"],atom_list,1)
  
#Save and overwrite active model:
def quicksave_active():
  import shutil
  mol_id = _active_molecule_or_status()
  if mol_id is None:
    return None
  filename=molecule_name(mol_id)
  filename_bak=filename+".bak"
  shutil.copy(filename,filename_bak)
  save_coordinates(mol_id,filename)
  info_dialog("Saved {filename} (molecule {mol_id}). \n If this was an accident, you can find a backup of the original file here: \n {filename_bak}".format(filename=filename,mol_id=mol_id,filename_bak=filename_bak))

#Refresh model from file (if changed externally for example)
def reload_model():
  mol_id = _active_molecule_or_status()
  if mol_id is None:
    return None
  filename=molecule_name(mol_id)
  clear_and_update_model_molecule_from_file(mol_id,filename)

#Copy changes from active chain to NCS equivalents.
def copy_ncs_chain_from_active():
  chain_context = _active_chain_context_or_status()
  if not chain_context:
    return None
  mol_id = chain_context["mol_id"]
  ch_id = chain_context["chain_id"]
  with _grouped_model_edit(mol_id):
    ncs_control_change_ncs_master_to_chain_id(mol_id, ch_id)
    copy_from_ncs_master_to_others(mol_id, ch_id)

def return_seq_as_string(mol_id,ch_id):
  seq_chars=[]
  for sn in range(0,chain_n_residues(ch_id,mol_id)):
    resno=seqnum_from_serial_number(mol_id,ch_id,sn)
    ins_code=insertion_code_from_serial_number(mol_id,ch_id,sn)
    residue_name_3_letter=residue_name(mol_id,ch_id,resno,ins_code)
    aa_code=three_letter_code2single_letter(residue_name_3_letter)
    if len(residue_name_3_letter)==1:
      aa_code=residue_name_3_letter
    elif residue_name_3_letter!="ALA" and aa_code=="A":
      aa_code="X"
    seq_chars.append(aa_code)
  return "".join(seq_chars)

def _compile_sequence_pattern(subseq):
  """Compile simple sequence search syntax such as X and (A/S)."""
  pattern=[]
  index=0
  while index < len(subseq):
    character=subseq[index]
    if character=="(":
      close_index=subseq.find(")", index+1)
      if close_index==-1:
        raise ValueError("Unmatched '(' in sequence pattern")
      option_text=subseq[index+1:close_index]
      if not option_text:
        raise ValueError("Empty bracketed option in sequence pattern")
      options=[piece.strip().upper() for piece in option_text.split("/") if piece.strip()]
      if not options:
        raise ValueError("Empty bracketed option in sequence pattern")
      allowed=set()
      for option in options:
        if len(option)!=1:
          raise ValueError("Bracketed options must be single-letter residues")
        allowed.add(option)
      pattern.append(("set", allowed))
      index=close_index+1
      continue
    if character=="X":
      pattern.append(("any", None))
    else:
      pattern.append(("set", set([character])))
    index=index+1
  if not pattern:
    raise ValueError("Empty sequence pattern")
  return pattern


def _sequence_pattern_matches_at(seq, start_index, compiled_pattern):
  if start_index+len(compiled_pattern) > len(seq):
    return False
  for offset in range(len(compiled_pattern)):
    mode, allowed=compiled_pattern[offset]
    seq_char=seq[start_index+offset]
    if mode=="any":
      continue
    if seq_char not in allowed:
      return False
  return True


def _sequence_match_context_label(ch_id, seq, sn_start, pattern_length, mol_id):
  sn_end=sn_start+pattern_length-1
  resno_start=seqnum_from_serial_number(mol_id,ch_id,sn_start)
  resno_end=seqnum_from_serial_number(mol_id,ch_id,sn_end)
  left_context_start=max(0, sn_start-3)
  right_context_end=min(len(seq), sn_start+pattern_length+3)
  left_context=seq[left_context_start:sn_start]
  match_context=seq[sn_start:sn_start+pattern_length]
  right_context=seq[sn_start+pattern_length:right_context_end]
  if left_context_start>0:
    left_context="..."+left_context
  if right_context_end<len(seq):
    right_context=right_context+"..."
  return "%s%d-%d: %s*%s*%s" %(ch_id,resno_start,resno_end,left_context,match_context,right_context)


def _sequence_match_list_entry(mol_id, ch_id, seq, sn_start, pattern_length):
  """Build a chooser entry for a sequence hit using the residue centre."""
  resno=seqnum_from_serial_number(mol_id,ch_id,sn_start)
  ins_code=insertion_code_from_serial_number(mol_id,ch_id,sn_start)
  residue_centre=residue_centre_py(mol_id,ch_id,resno,ins_code)
  if not (isinstance(residue_centre, (list, tuple)) and len(residue_centre)==3):
    return None
  return [
    _sequence_match_context_label(ch_id,seq,sn_start,pattern_length,mol_id),
    residue_centre[0],
    residue_centre[1],
    residue_centre[2],
    {
      "type": "polymer_residue",
      "mol_id": mol_id,
      "chain_id": ch_id,
      "resno": resno,
      "ins_code": ins_code or "",
      "serial_number": sn_start,
    },
  ]


def find_sequence_in_current_chain(subseq):
  subseq=subseq.upper().strip()
  residue=_active_residue_or_status()
  if not residue:
    return None
  mol_id=residue[0]
  ch_id=residue[1]
  seq=return_seq_as_string(mol_id,ch_id)
  try:
    compiled_pattern=_compile_sequence_pattern(subseq)
  except ValueError as error:
    info_dialog(str(error))
    return None
  sn_list=[]
  interesting_list=[]
  pattern_length=len(compiled_pattern)
  index=0
  while index <= len(seq)-pattern_length:
    if _sequence_pattern_matches_at(seq, index, compiled_pattern):
      sn_list.append(index)
      index=index+pattern_length
    else:
      index=index+1
  if len(sn_list)==0:
    info_dialog("Sequence not found!")
    return None
  for sn in sn_list:
    list_entry=_sequence_match_list_entry(mol_id,ch_id,seq,sn,pattern_length)
    if not list_entry:
      continue
    interesting_list.append(list_entry)
  if not interesting_list:
    if len(sn_list)==1:
      info_dialog("Found the sequence match, but could not locate a CA/P atom to center on.")
    else:
      info_dialog("Found matching sequence positions, but could not locate jump atoms for them.")
    return None
  result_label = "Match to entered sequence" if len(interesting_list) == 1 else "Matches to entered sequence"
  navigable_interesting_things_gui(result_label, interesting_list)
  return None

# ---------------------------------------------------------------------------
# Legacy Gtk archive
# ---------------------------------------------------------------------------
# The primary Coot 1.x startup/menu path lives above. The large blocks below
# preserve older Gtk-heavy helpers so we can either revive them selectively or
# refer back to the original implementations when porting features.
#
# Historical dialog-backed helper implementations promoted back into normal
# Python. The remaining legacy notes below are kept only as brief reference.
#
# Archived startup wiring that is currently non-functional in this build:
#
# add_module_cryo_em_gui()
# add_module_refine()
# #add_module_restraints()
# if GUI_PYTHON_AVAILABLE:
#   add_module_carbohydrate_gui()
#
# coot_toolbar_button("Measure distance", "do_distance_define()", icon_name="geom.svg")
# coot_toolbar_button("Sym?", "set_show_symmetry_master(not get_show_symmetry())", icon_name="cell+symm.svg")
# _register_key_binding_if(
#   DIALOG_KEYBINDINGS_AVAILABLE,
#   "Place atom at pointer",
#   "P",
#   lambda: place_atom_at_pointer(),
# )
# _register_key_binding_if(
#   DIALOG_KEYBINDINGS_AVAILABLE,
#   "Set map contour in sigma",
#   "L",
#   lambda: set_map_level_quickly(),
# )
# coot_toolbar_button("Sequence context", "sequence_context()", icon_name="")
# coot_toolbar_button("Accept RSR", "accept_regularizement()", icon_name="")
#
# Historical dialog-backed helper implementations promoted back into normal
# Python code.
def set_map_level_quickly():
  map_id = _scrollable_map_or_status()
  if map_id is None:
    return None
  if map_is_displayed(map_id)==0:
    info_dialog("You need a (scrollable, displayed) map!")
    return None
  current_map_level="{0:.2f}".format(get_contour_level_in_sigma(map_id))
  def set_map_level_quickly_from_entry(X):
    try:
      map_level=float(str(X).strip())
    except ValueError:
      info_dialog("Has to be a number!")
      return 0
    set_contour_level_in_sigma(map_id, map_level)
    _post_map_level_status(map_id, "set to")
    return 1
  generic_single_entry("New map level in sigma/RMS?",current_map_level,"Set map level",set_map_level_quickly_from_entry)

def _grow_secondary_structure_from_click(prompt_label, button_label, phi, psi, set_b_factor_range=True):
  """Grow a secondary-structure segment from the clicked terminus."""
  def grow_post_click(res1):
    mol_id=_click_spec_imol(res1)
    ch_id=_click_spec_chain_id(res1)
    res_no_start=_click_spec_res_no(res1)
    if mol_id == -1 or ch_id is False or res_no_start is False:
      info_dialog(f"{button_label} requires a valid clicked residue.")
      return 0

    def grow_from_entry(n):
      count = _parse_int_entry(n, f"{button_label} requires a positive integer.", minimum_value=1)
      if count is None:
        return 0
      res_no=res_no_start
      res_no_0=res_no
      for _i in range(1, count+1):
        add_terminal_residue_using_phi_psi(mol_id, ch_id, res_no, "auto", phi, psi)
        sort_residues(mol_id)
        if res_no == (first_residue_in_seg(mol_id, ch_id, res_no)+1):
          res_no=res_no-1
        elif res_no == (last_residue_in_seg(mol_id, ch_id, res_no)-1):
          res_no=res_no+1
      if set_b_factor_range:
        set_b_factor_residue_range(mol_id, ch_id, res_no_0, res_no, default_new_atoms_b_factor())
      return 1

    generic_single_entry(prompt_label, "10", button_label, grow_from_entry)
  _script_user_defined_click(1, grow_post_click)


#Grow helix from selected terminus
def grow_helix():
  return _grow_secondary_structure_from_click(
    "How many residues for helix?",
    "Grow helix",
    -57.82,
    -47,
    set_b_factor_range=True,
  )
  
#Grow strand from selected terminus
def grow_strand():
  return _grow_secondary_structure_from_click(
    "How many residues for strand?",
    "Grow strand",
    -139,
    135,
    set_b_factor_range=True,
  )
  
#Grow para strand from selected terminus
def grow_parallel_strand():
  return _grow_secondary_structure_from_click(
    "How many residues for parallel strand?",
    "Grow parallel strand",
    -119,
    113,
    set_b_factor_range=False,
  )

#Grow 3-10 helix from selected terminus
def grow_helix_3_10():
  return _grow_secondary_structure_from_click(
    "How many residues for helix?",
    "Grow 3-10 helix",
    -49,
    -26,
    set_b_factor_range=False,
  )

#Renumbers the active chain (from active_residue()). User enters new starting residue number.
def renumber_by_first_res():
  chain_context = _active_chain_context_or_status()
  if not chain_context:
    return None
  mol_id=chain_context["mol_id"]
  ch_id=chain_context["chain_id"]
  first_res=first_residue(mol_id,ch_id)
  last_res=last_residue(mol_id,ch_id)
  def renum_chain(new_num):
    new_num = _parse_int_entry(new_num, "Renumbering requires an integer residue number.")
    if new_num is None:
      return 0
    offset=new_num-first_res
    if not _renumber_residue_range_or_dialog(mol_id, ch_id, first_res, last_res, offset):
      return 0
    _refresh_extra_restraints_display(mol_id)
    return 1
  generic_single_entry("New number for first residue in this chain?",
  str(first_res),"Renumber",renum_chain)
  
#Renumbers the active chain (from active_residue()). User enters new last residue number.
def renumber_by_last_res():
  chain_context = _active_chain_context_or_status()
  if not chain_context:
    return None
  mol_id=chain_context["mol_id"]
  ch_id=chain_context["chain_id"]
  first_res=first_residue(mol_id,ch_id)
  last_res=last_residue(mol_id,ch_id)
  def renum_chain(new_num):
    new_num = _parse_int_entry(new_num, "Renumbering requires an integer residue number.")
    if new_num is None:
      return 0
    offset=new_num-last_res
    if not _renumber_residue_range_or_dialog(mol_id, ch_id, first_res, last_res, offset):
      return 0
    _refresh_extra_restraints_display(mol_id)
    return 1
  generic_single_entry("New number for last residue in this chain?",
  str(last_res),"Renumber",renum_chain)
  
#Renumbers the active chain (from active_residue()). User enters new active residue number.
def renumber_by_active_res():
  chain_context = _active_chain_context_or_status()
  if not chain_context:
    return None
  mol_id=chain_context["mol_id"]
  ch_id=chain_context["chain_id"]
  current_num=chain_context["resno"]
  first_res=first_residue(mol_id,ch_id)
  last_res=last_residue(mol_id,ch_id)
  def renum_chain(new_num):
    new_num = _parse_int_entry(new_num, "Renumbering requires an integer residue number.")
    if new_num is None:
      return 0
    offset=new_num-current_num
    if not _renumber_residue_range_or_dialog(mol_id, ch_id, first_res, last_res, offset):
      return 0
    _refresh_extra_restraints_display(mol_id)
    return 1
  generic_single_entry("New number for this residue?",
  str(current_num),"Renumber",renum_chain)

#Renumber from N-term to current residue
def renumber_n_term_segment():
  chain_context = _active_chain_context_or_status()
  if not chain_context:
    return None
  mol_id=chain_context["mol_id"]
  ch_id=chain_context["chain_id"]
  resn=chain_context["resno"]
  first_res=seqnum_from_serial_number(mol_id,"%s"%(ch_id),0)
  def renumber_n_term_segment_entry(new_resn):
    new_resn = _parse_int_entry(new_resn, "Renumbering requires an integer residue number.")
    if new_resn is None:
      return 0
    offset=new_resn-resn
    if not _renumber_residue_range_or_dialog(mol_id, ch_id, first_res, resn, offset):
      return 0
    _refresh_extra_restraints_display(mol_id)
    return 1
  generic_single_entry("New residue number?",
  str(resn),"Renumber",renumber_n_term_segment_entry)
  
#Renumber from current residue to C-term
def renumber_c_term_segment():
  chain_context = _active_chain_context_or_status()
  if not chain_context:
    return None
  mol_id=chain_context["mol_id"]
  ch_id=chain_context["chain_id"]
  resn=chain_context["resno"]
  last_res=last_residue(mol_id,ch_id)
  def renumber_c_term_segment_entry(new_resn):
    new_resn = _parse_int_entry(new_resn, "Renumbering requires an integer residue number.")
    if new_resn is None:
      return 0
    offset=new_resn-resn
    if not _renumber_residue_range_or_dialog(mol_id, ch_id, resn, last_res, offset):
      return 0
    _refresh_extra_restraints_display(mol_id)
    return 1
  generic_single_entry("New residue number?",
  str(resn),"Renumber",renumber_c_term_segment_entry)
  
#Reload map with different hi-res limit
def change_hires_limit():
  map_id = _scrollable_map_or_status()
  if map_id is None:
    info_dialog("You need a map!")
  else:
    def change_hires_by_entry(new_res):
      try:
        new_res=float(str(new_res).strip())
      except Exception:
        info_dialog("High-resolution limit must be a number.")
        return 0
      mtz_file=map_parameters(map_id)[0]
      F_col=map_parameters(map_id)[1]
      PHI_col=map_parameters(map_id)[2]
      make_and_draw_map_with_reso_with_refmac_params(mtz_file,F_col,PHI_col,"",0,0,0,"Fobs:None-specified",
      "SigF:None-specified","RFree:None-specified",0,0,1,1000.0,new_res)
      close_molecule(map_id)
      return 1
    generic_single_entry("New high-res limit for map?",
    "5.0","Change high resolution limit for active map",change_hires_by_entry)
    


import hashlib
import json
import os
import shutil
import traceback
from functools import lru_cache


COMMON_MONOMER_FAVORITES_FILENAME = "coot_trimmings_favorites.json"
COMMON_MONOMER_FAVORITE_CIF_PREFIX = "coot_trimmings_favorite_"
COMMON_MONOMER_FAVORITES_MENU = None

COMMON_MONOMER_MENU = [
  ("Buffers", [
    ("Bicine (BCN)", "BCN"),
    ("Bis-Tris (BTB)", "BTB"),
    ("Cacodylate (CAC)", "CAC"),
    ("HEPES (EPE)", "EPE"),
    ("MES (MES)", "MES"),
    ("TAPS (T3A)", "T3A"),
    ("Triethanolamine (TAM)", "TAM"),
    ("Tris (TRS)", "TRS"),
  ]),
  ("Detergents / lipids", [
    ("Detergents", [
      ("Alpha-DDM / dodecyl-alpha-maltoside (LMU)", "LMU"),
      ("BGL / beta-2-octylglucoside (BGL)", "BGL"),
      ("BNG / beta-nonylglucoside (BNG)", "BNG"),
      ("BOG / beta-octylglucoside (BOG)", "BOG"),
      ("C10E6 (C10)", "C10"),
      ("C8E4 (C8E)", "C8E"),
      ("CHAPS (CPS)", "CPS"),
      ("CHAPSO (1N7)", "1N7"),
      ("Cyclohexyl-hexyl-beta-D-maltoside (MA4)", "MA4"),
      ("CYMAL-4 (CVM)", "CVM"),
      ("CYMAL-5 (CM5)", "CM5"),
      ("DDM (dodecyl maltoside) (LMT)", "LMT"),
      ("Digitonin (AJP)", "AJP"),
      ("DM (decyl maltoside) (DMU)", "DMU"),
      ("LDAO (LDA)", "LDA"),
      ("LMNG (LMN)", "LMN"),
      ("Octyl beta-D-galactopyranoside (HSH)", "HSH"),
      ("OGNG (37X)", "37X"),
      ("Seleno-DDM / dodecyl-beta-selenomaltoside (LSM)", "LSM"),
      ("Undecyl maltoside (UMQ)", "UMQ"),
    ]),
    ("Lipids", [
      ("Bile salts", [
        ("Cholic acid (CHD)", "CHD"),
        ("Deoxycholic acid (DXC)", "DXC"),
        ("Taurodeoxycholate (6SB)", "6SB"),
      ]),
      ("Dicaproyl phosphatidylserine (PSF)", "PSF"),
      ("Dioleoyl phosphatidylcholine (PCW)", "PCW"),
      ("Dipalmitoyl phosphatidic acid (PX6)", "PX6"),
      ("Dipalmitoyl phosphatidylethanolamine (PEF)", "PEF"),
      ("Dipalmitoyl phosphatidylglycerol (LHG)", "LHG"),
      ("Diundecyl phosphatidylcholine (PLC)", "PLC"),
      ("PIPs / inositides", [
        ("Dihexadecanoyl PI(3,4,5)P3 (PIZ)", "PIZ"),
        ("Dihexadecanoyl PI(4,5)P2 (PIK)", "PIK"),
        ("Dioctanoyl PI(3,4)P2 (52N)", "52N"),
        ("Dioctanoyl PI(3,4,5)P3 (IP9)", "IP9"),
        ("Phosphatidylinositol (T7X)", "T7X"),
      ]),
      ("Sphingolipids", [
        ("Sphingomyelin (FO4)", "FO4"),
        ("Sphingosine (SPH)", "SPH"),
        ("Sphingosine-1-phosphate (S1P)", "S1P"),
        ("Sulfatide (SLF)", "SLF"),
        ("Sulfogalactocerebroside (SFT)", "SFT"),
      ]),
    ]),
    ("Sterols", [
      ("Cholesterol (CLR)", "CLR"),
      ("Cholesteryl hemisuccinate (Y01)", "Y01"),
    ]),
  ]),
  ("Ions / metals", [
    ("Cations", [
      ("Ammonium (NH4)", "NH4"),
      ("Barium (BA)", "BA"),
      ("Calcium (CA)", "CA"),
      ("Cesium (CS)", "CS"),
      ("Lithium (LI)", "LI"),
      ("Magnesium (MG)", "MG"),
      ("Potassium (K)", "K"),
      ("Rubidium (RB)", "RB"),
      ("Sodium (NA)", "NA"),
      ("Strontium (SR)", "SR"),
    ]),
    ("Halides", [
      ("Bromide (BR)", "BR"),
      ("Chloride (CL)", "CL"),
      ("Fluoride (F)", "F"),
      ("Iodide (IOD)", "IOD"),
    ]),
    ("Heavy atoms", [
      ("Europium(III) ion (EU3)", "EU3"),
      ("Gadolinium ion (GD3)", "GD3"),
      ("Gold ion (AU)", "AU"),
      ("Holmium(III) ion (HO3)", "HO3"),
      ("Iridium ion (IR)", "IR"),
      ("Krypton (KR)", "KR"),
      ("Lanthanum(III) ion (LA)", "LA"),
      ("Lead(II) ion (PB)", "PB"),
      ("Lutetium(III) ion (LU)", "LU"),
      ("Mercury(II) ion (HG)", "HG"),
      ("Osmium ion (OS)", "OS"),
      ("Platinum(II) ion (PT)", "PT"),
      ("Samarium(III) ion (SM)", "SM"),
      ("Terbium(III) ion (TB)", "TB"),
      ("Uranium (U1)", "U1"),
      ("Uranyl(VI) ion (IUM)", "IUM"),
      ("Xenon (XE)", "XE"),
      ("Ytterbium(III) ion (YB)", "YB"),
    ]),
    ("Metal clusters", [
      ("Fe-S cluster (35L)", "35L"),
      ("Fe-S cluster (FS1)", "FS1"),
      ("Fe-S cluster / cubane (SF4)", "SF4"),
      ("FeMo cofactor-like cluster (ICS)", "ICS"),
      ("Ni-Fe cluster (82N)", "82N"),
      ("Tetranuclear copper-sulfide cluster (CUZ)", "CUZ"),
    ]),
    ("Other anions", [
      ("Azide (AZI)", "AZI"),
      ("Biselenite (BSY)", "BSY"),
      ("Carbonate (CO3)", "CO3"),
      ("Nitrate (NO3)", "NO3"),
      ("Perchlorate (LCP)", "LCP"),
      ("Phosphite (PO3)", "PO3"),
      ("Sulfite (SO3)", "SO3"),
      ("Thiocyanate (SCN)", "SCN"),
    ]),
    ("Oxyanions", [
      ("Arsenate (ART)", "ART"),
      ("Molybdate (MOO)", "MOO"),
      ("Phosphate (PO4)", "PO4"),
      ("Selenate (SE4)", "SE4"),
      ("Sulfate (SO4)", "SO4"),
      ("Tungstate (WO4)", "WO4"),
      ("Vanadate (VO4)", "VO4"),
    ]),
    ("Transition / heavy metals", [
      ("Barium (BA)", "BA"),
      ("Cadmium (CD)", "CD"),
      ("Cobalt (CO)", "CO"),
      ("Copper (CU)", "CU"),
      ("Gold ion (AU)", "AU"),
      ("Iron (FE)", "FE"),
      ("Manganese (MN)", "MN"),
      ("Mercury(II) ion (HG)", "HG"),
      ("Nickel (NI)", "NI"),
      ("Platinum(II) ion (PT)", "PT"),
      ("Strontium (SR)", "SR"),
      ("Zinc (ZN)", "ZN"),
    ]),
  ]),
  ("Ligands", [
    ("Bases / nucleosides", [
      ("Adenine (ADE)", "ADE"),
      ("Cytosine (CYT)", "CYT"),
      ("Guanosine (GMP)", "GMP"),
      ("Uracil (URA)", "URA"),
      ("Uridine (URI)", "URI"),
    ]),
    ("Cofactors", [
      ("Acidic NAD+ form (NAJ)", "NAJ"),
      ("Adenosine 5'-pentaphosphate (5FA)", "5FA"),
      ("ADP-ribose (APR)", "APR"),
      ("Coenzyme A (COA)", "COA"),
      ("Cyclic ADP-ribose (CXR)", "CXR"),
      ("Diadenosine pentaphosphate (AP5)", "AP5"),
      ("Diadenosine tetraphosphate (B4P)", "B4P"),
      ("Diadenosine triphosphate (BA3)", "BA3"),
      ("Diguanosine pentaphosphate (GP5)", "GP5"),
      ("Diguanosine triphosphate (GP3)", "GP3"),
      ("Flavin adenine dinucleotide (FAD)", "FAD"),
      ("Flavin mononucleotide (FMN)", "FMN"),
      ("Nicotinamide adenine dinucleotide (NAD)", "NAD"),
      ("Nicotinamide adenine dinucleotide phosphate (NAP)", "NAP"),
      ("Nicotinic acid adenine dinucleotide phosphate (DN4)", "DN4"),
      ("ppGpp / guanosine 5',3'-tetraphosphate (G4P)", "G4P"),
      ("Reduced nicotinamide adenine dinucleotide (NAI)", "NAI"),
      ("Reduced nicotinamide adenine dinucleotide phosphate (NDP)", "NDP"),
      ("S-adenosylhomocysteine (SAH)", "SAH"),
      ("S-adenosylmethionine (SAM)", "SAM"),
    ]),
    ("Free amino acids", [
      ("Alanine (ALA)", "ALA"),
      ("Arginine (ARG)", "ARG"),
      ("Asparagine (ASN)", "ASN"),
      ("Aspartic acid (ASP)", "ASP"),
      ("Cysteine (CYS)", "CYS"),
      ("Glutamic acid (GLU)", "GLU"),
      ("Glutamine (GLN)", "GLN"),
      ("Glycine (GLY)", "GLY"),
      ("Histidine (HIS)", "HIS"),
      ("Isoleucine (ILE)", "ILE"),
      ("Leucine (LEU)", "LEU"),
      ("Lysine (LYS)", "LYS"),
      ("Methionine (MET)", "MET"),
      ("Phenylalanine (PHE)", "PHE"),
      ("Proline (PRO)", "PRO"),
      ("Serine (SER)", "SER"),
      ("Threonine (THR)", "THR"),
      ("Tryptophan (TRP)", "TRP"),
      ("Tyrosine (TYR)", "TYR"),
      ("Valine (VAL)", "VAL"),
    ]),
    ("Free sugars", [
      ("Alpha-D-galactose (GLA)", "GLA"),
      ("Alpha-D-glucose (GLC)", "GLC"),
      ("Alpha-D-ribofuranose (RIB)", "RIB"),
      ("Alpha-L-arabinose (ARA)", "ARA"),
      ("Beta-D-galactose (GAL)", "GAL"),
      ("Beta-D-glucose (BGC)", "BGC"),
      ("Fructose (FRU)", "FRU"),
      ("Fucose (FUC)", "FUC"),
      ("Mannose (MAN)", "MAN"),
      ("Xylose (XYS)", "XYS"),
    ]),
    ("Hemes / porphyrins", [
      ("Ferrous tetravinylporphine complex (HEV)", "HEV"),
      ("Heme A (HEA)", "HEA"),
      ("Heme B (HEM)", "HEM"),
      ("Heme C (HEC)", "HEC"),
      ("Heme D (DHE)", "DHE"),
      ("Heme D hydroxychlorin spirolactone (HDD)", "HDD"),
      ("Heme O (HEO)", "HEO"),
      ("Heme-AS (HAS)", "HAS"),
    ]),
    ("Immunosuppressants / macrocycles", [
      ("Ascomycin analog (FK5)", "FK5"),
      ("Rapamycin (RAP)", "RAP"),
    ]),
    ("Non-hydrolysable nucleotide analogs", [
      ("AMP-CPP (APC)", "APC"),
      ("AMP-PCP (ACP)", "ACP"),
      ("AMP-PNP / AMPPNP (ANP)", "ANP"),
      ("GMP-PCP (GCP)", "GCP"),
      ("GMP-PNP / GppNHp (GNP)", "GNP"),
      ("GTPgammaS (GSP)", "GSP"),
    ]),
    ("Nucleotides", [
      ("Adenosine diphosphate (ADP)", "ADP"),
      ("Adenosine monophosphate (AMP)", "AMP"),
      ("Adenosine triphosphate (ATP)", "ATP"),
      ("Cytidine diphosphate (CDP)", "CDP"),
      ("Cytidine monophosphate (C5P)", "C5P"),
      ("Cytidine triphosphate (CTP)", "CTP"),
      ("Guanosine diphosphate (GDP)", "GDP"),
      ("Guanosine monophosphate (5GP)", "5GP"),
      ("Guanosine triphosphate (GTP)", "GTP"),
      ("Uridine diphosphate (UDP)", "UDP"),
      ("Uridine monophosphate (U)", "U"),
      ("Uridine triphosphate (UTP)", "UTP"),
    ]),
    ("Phosphosugars / inositol phosphates", [
      ("Fructose-1,6-bisphosphate (FBP)", "FBP"),
      ("Fructose-6-phosphate (F6P)", "F6P"),
      ("Glucosamine-6-phosphate (GLP)", "GLP"),
      ("Glucose-1-phosphate (G1P)", "G1P"),
      ("Glucose-6-phosphate (G6P)", "G6P"),
      ("Inositol hexakisphosphate / phytate (IHP)", "IHP"),
      ("Inositol pentakisphosphate (I5P)", "I5P"),
      ("Inositol tetrakisphosphate (4IP)", "4IP"),
      ("Inositol-1,4,5-trisphosphate (I3P)", "I3P"),
      ("Inositol-1,4-bisphosphate (2IP)", "2IP"),
      ("Inositol-1-phosphate (IPD)", "IPD"),
      ("Inositol-2,4,5-trisphosphate (I2P)", "I2P"),
      ("Inositol-4,5-bisphosphate (IP2)", "IP2"),
      ("Inositol-4-phosphate (I4D)", "I4D"),
      ("Mannose-1-phosphate (M1P)", "M1P"),
      ("Mannose-6-phosphate (M6P)", "M6P"),
      ("N-acetylglucosamine-6-phosphate (4QY)", "4QY"),
      ("Ribose-5-phosphate (R5P)", "R5P"),
    ]),
    ("Polyamines", [
      ("Cadaverine (N2P)", "N2P"),
      ("Putrescine (PUT)", "PUT"),
      ("Spermidine (SPD)", "SPD"),
      ("Spermine (SPM)", "SPM"),
    ]),
    ("Protease inhibitors", [
      ("AEBSF (AES)", "AES"),
      ("Benzamidine (BEN)", "BEN"),
      ("E-64 (E64)", "E64"),
      ("PMSF (PMF)", "PMF"),
      ("TLCK (TCK)", "TCK"),
    ]),
    ("Retinoids", [
      ("9-cis-retinoic acid (9CR)", "9CR"),
      ("All-trans retinoic acid (REA)", "REA"),
      ("Retinal (RET)", "RET"),
      ("Retinol (RTL)", "RTL"),
    ]),
    ("Tetrapyrroles / chlorophylls", [
      ("Bacteriochlorophyll A (BCL)", "BCL"),
      ("Bacteriochlorophyll B (BCB)", "BCB"),
      ("Biliverdin (EL5)", "EL5"),
      ("Chlorophyll A (CLA)", "CLA"),
      ("Chlorophyll B (CHL)", "CHL"),
      ("Chlorophyll D (CL7)", "CL7"),
      ("Chlorophyll F (F6C)", "F6C"),
      ("Pheophytin A (PHO)", "PHO"),
      ("Phycocyanobilin (CYC)", "CYC"),
      ("Protoporphyrin IX (PP9)", "PP9"),
    ]),
    ("Vitamin / coenzyme cofactors", [
      ("5,10-Methenyltetrahydrofolate (GUE)", "GUE"),
      ("Biotin (BTN)", "BTN"),
      ("Pyridoxal 5'-phosphate (PLP)", "PLP"),
      ("Tetrahydrofolate (THG)", "THG"),
      ("Thiamine diphosphate (TPP)", "TPP"),
    ]),
  ]),
  ("Solvents / additives", [
    ("1,4-Butanediol (BU1)", "BU1"),
    ("1-Butanol (1BO)", "1BO"),
    ("Acetate (ACT)", "ACT"),
    ("Acetonitrile (CCN)", "CCN"),
    ("Ammonium (NH4)", "NH4"),
    ("Azide (AZI)", "AZI"),
    ("Carbonate (CO3)", "CO3"),
    ("Chelators", [
      ("EDTA (EDT)", "EDT"),
      ("Nitrilotriacetic acid (NTA)", "NTA"),
    ]),
    ("DMSO (DMS)", "DMS"),
    ("Ethanol (EOH)", "EOH"),
    ("Ethylene glycol (EDO)", "EDO"),
    ("Formic acid (FMT)", "FMT"),
    ("Glycerol (GOL)", "GOL"),
    ("Isopropanol (IPA)", "IPA"),
    ("Methanol (MOH)", "MOH"),
    ("MPD (MPD)", "MPD"),
    ("n-Propanol (POL)", "POL"),
    ("Nitrate (NO3)", "NO3"),
    ("Organic acids", [
      ("Citrate anion (FLC)", "FLC"),
      ("Citric acid (CIT)", "CIT"),
      ("Isocitric acid (ICT)", "ICT"),
      ("Malate (MLT)", "MLT"),
      ("Malonate (MLI)", "MLI"),
      ("Oxalate (OXL)", "OXL"),
      ("Oxaloacetate (OAA)", "OAA"),
    ]),
    ("Organic cations", [
      ("Tetrabutylammonium (TBA)", "TBA"),
      ("Tetramethylammonium (TMA)", "TMA"),
      ("Triethylammonium (TEA)", "TEA"),
    ]),
    ("PEG (PEG)", "PEG"),
    ("Propylene glycol / R-1,2-propanediol (PGR)", "PGR"),
    ("Propylene glycol / S-1,2-propanediol (PGO)", "PGO"),
    ("Reducing agents", [
      ("Beta-mercaptoethanol (BME)", "BME"),
      ("DTT / 2,3-dihydroxy-1,4-dithiobutane (DTT)", "DTT"),
      ("TCEP / tris(2-carboxyethyl)phosphine (TCE)", "TCE"),
    ]),
    ("Tartaric acid", [
      ("D(-)-Tartaric acid (TAR)", "TAR"),
      ("L(+)-Tartaric acid (TLA)", "TLA"),
    ]),
    ("Taurine (TAU)", "TAU"),
    ("Tetraethylene glycol (PG4)", "PG4"),
    ("Thiocyanate (SCN)", "SCN"),
    ("TMAO / trimethylamine oxide (TMO)", "TMO"),
    ("Urea (URE)", "URE"),
  ]),
]

# Some "monomers" are really better handled as typed atoms or built-in simple
# groups. Those are placed at the pointer instead of going through get_monomer().
COMMON_MONOMER_POINTER_TYPES = {
  "H": "H",
  "HE": "HE",
  "LI": "LI",
  "BE": "BE",
  "B": "B",
  "C": "C",
  "N": "N",
  "O": "O",
  "F": "F",
  "CL": "CL",
  "BR": "BR",
  "IOD": "IOD",
  "ZN": "ZN",
  "MN": "MN",
  "FE": "FE",
  "CO": "CO",
  "CU": "CU",
  "NI": "NI",
  "CD": "CD",
  "HG": "HG",
  "PB": "PB",
  "AU": "AU",
  "PT": "PT",
  "OS": "OS",
  "IR": "IR",
  "GD3": "GD3",
  "EU3": "EU3",
  "TB": "TB",
  "HO3": "HO3",
  "LA": "LA",
  "LU": "LU",
  "SM": "SM",
  "YB": "YB",
  "XE": "XE",
  "KR": "KR",
  "U1": "U1",
  "PO4": "PO4",
  "SO4": "SO4",
}



def _gio_module():
  gio_module = globals().get("Gio")
  coot_gui_module = _coot_gui_module()
  if gio_module is None and coot_gui_module is not None and hasattr(coot_gui_module, "Gio"):
    gio_module = coot_gui_module.Gio
    globals()["Gio"] = gio_module
  return gio_module


def _common_monomer_favorites_path():
  return os.path.join(_coot_trimmings_dir(), COMMON_MONOMER_FAVORITES_FILENAME)


RECENT_FILES_LIMIT = 20
RECENT_FILES_FILENAME = "coot_trimmings_recent_files.json"


def _recent_files_path():
  return os.path.join(_coot_trimmings_dir(), RECENT_FILES_FILENAME)


def _recent_file_normalize(entry, category):
  if not isinstance(entry, dict) or not isinstance(entry.get("path"), str) or not entry["path"].strip():
    return None
  path = os.path.realpath(os.path.abspath(os.path.expanduser(entry["path"])))
  kind = entry.get("kind")
  if kind not in (("model",) if category == "models" else ("map", "mtz")):
    return None
  result = {"path": path, "kind": kind}
  if kind == "map":
    result["difference"] = bool(entry.get("difference", False))
  if kind == "mtz":
    maps = entry.get("maps")
    if not isinstance(maps, list):
      return None
    result["maps"] = []
    for row in maps[:32]:
      if not isinstance(row, dict) or any(not isinstance(row.get(key), str) or not row[key].strip()
                                         for key in ("f", "phi")):
        continue
      if not isinstance(row.get("weight", ""), str):
        continue
      normalized = {"f": row["f"], "phi": row["phi"], "weight": row.get("weight", ""),
                    "use_weights": bool(row.get("use_weights", False)),
                    "difference": bool(row.get("difference", False))}
      if normalized["use_weights"] and not normalized["weight"].strip():
        continue
      if normalized not in result["maps"]:
        result["maps"].append(normalized)
    if not result["maps"]:
      return None
  return result


def _recent_files_state():
  state = getattr(coot, "_trimmings_recent_files", None)
  if state is not None:
    return state
  state = {"data": {"models": [], "maps": []}, "seen": {}, "menus": {}, "timer": None,
           "dirty": False, "write_blocked": False, "error": None}
  try:
    with open(_recent_files_path(), "r", encoding="utf-8") as handle:
      raw = json.load(handle)
    if not isinstance(raw, dict) or raw.get("version") != 1:
      raise ValueError("Unsupported recent-files format")
    for category in ("models", "maps"):
      rows = raw.get(category, [])
      if not isinstance(rows, list):
        raise ValueError("Invalid recent-files list")
      seen = set()
      for row in rows:
        entry = _recent_file_normalize(row, category)
        if entry is not None and entry["path"] not in seen:
          state["data"][category].append(entry)
          seen.add(entry["path"])
        if len(state["data"][category]) == RECENT_FILES_LIMIT:
          break
  except FileNotFoundError:
    pass
  except (OSError, ValueError, TypeError) as error:
    state["write_blocked"] = True
    print("Recent files: could not load history; preserving the existing file:", error)
  coot._trimmings_recent_files = state
  return state


def _save_recent_files():
  state = _recent_files_state()
  if not state["dirty"] or state["write_blocked"]:
    return False
  temporary = None
  try:
    path = _recent_files_path()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=".coot_recent_", suffix=".json", dir=os.path.dirname(path))
    with os.fdopen(fd, "w", encoding="utf-8") as handle:
      json.dump({"version": 1, **state["data"]}, handle, indent=2, ensure_ascii=True)
      handle.write("\n")
    os.replace(temporary, path)
    temporary = None
    state["dirty"], state["error"] = False, None
    return True
  except OSError as error:
    if state["error"] != str(error):
      print("Recent files: history is available for this session but could not be saved:", error)
      state["error"] = str(error)
    return False
  finally:
    if temporary is not None:
      try:
        os.remove(temporary)
      except OSError:
        pass


def _record_recent_file(entry, category):
  entry = _recent_file_normalize(entry, category)
  if entry is None:
    return False
  state = _recent_files_state()
  rows = state["data"][category]
  previous = next((row for row in rows if row["path"] == entry["path"]), None)
  if previous and entry["kind"] == previous["kind"] == "mtz":
    # One MTZ entry restores all coefficient sets observed for that file.
    entry["maps"] += [row for row in previous["maps"] if row not in entry["maps"]]
    entry["maps"] = entry["maps"][:32]
  updated = [entry] + [row for row in rows if row["path"] != entry["path"]]
  updated = updated[:RECENT_FILES_LIMIT]
  if updated == rows:
    return False
  state["data"][category] = updated
  state["dirty"] = True
  return True


def _recent_molecule_metadata(imol):
  if coot.is_valid_model_molecule(imol):
    return ("models", coot.molecule_name(imol), None, False)
  if coot.is_valid_map_molecule(imol):
    getter = getattr(coot, "map_parameters_py", None)
    params = getter(imol) if callable(getter) else None
    if isinstance(params, (list, tuple)):
      params = tuple(params)
    return ("maps", coot.molecule_name(imol), params, bool(coot.map_is_difference_map(imol)))
  return None


def _recent_entry_from_metadata(metadata):
  category, name, params, difference = metadata
  if category == "maps" and isinstance(params, (list, tuple)) and len(params) == 5 and params[0]:
    entry = _recent_file_normalize({"path": params[0], "kind": "mtz", "maps": [
      {"f": params[1], "phi": params[2], "weight": params[3], "use_weights": params[4],
       "difference": difference}]}, category)
  else:
    entry = _recent_file_normalize({"path": name, "kind": "model" if category == "models" else "map",
                                    "difference": difference}, category)
  # Derived maps, fragments and monomers normally have names, not file paths.
  if entry is None or not os.path.isfile(entry["path"]):
    return None
  if category == "maps" and entry["kind"] == "map" and entry["path"].lower().endswith(".mtz"):
    return None  # Never guess map columns when the build cannot report them.
  return entry


def _poll_recent_files():
  state = _recent_files_state()
  changed = False
  current = {}
  try:
    for imol in range(coot.graphics_n_molecules()):
      try:
        metadata = _recent_molecule_metadata(imol)
      except Exception as error:
        message = str(error)
        if state.get("poll_error") != message:
          print("Recent files: could not inspect molecule #%d:" % imol, error)
          state["poll_error"] = message
        continue
      if metadata is None:
        continue
      current[imol] = metadata
      if metadata == state["seen"].get(imol):
        continue
      entry = _recent_entry_from_metadata(metadata)
      if entry is not None:
        changed = _record_recent_file(entry, metadata[0]) or changed
    state["seen"] = current
    if changed or state["dirty"]:
      _save_recent_files()
    if changed:
      _refresh_recent_files_menus()
  except Exception as error:
    print("Recent files: metadata check failed:", error)
  return True


def _recent_short_directory(directory, tail=2):
  directory = os.path.normpath(directory)
  home = os.path.normpath(os.path.expanduser("~"))
  if directory == home:
    return "~"
  if directory.startswith(home + os.sep):
    prefix, parts = "~/", directory[len(home)+1:].split(os.sep)
    head = 1
  else:
    prefix, parts = os.sep, directory.strip(os.sep).split(os.sep)
    head = 2 if parts[0] == "Volumes" else 1
  if len(parts) > head+tail:
    parts = parts[:head] + ["..."] + parts[-tail:]
  return prefix + os.sep.join(parts)


def _recent_file_label(entry, rows):
  name = os.path.basename(entry["path"])
  directory = os.path.dirname(entry["path"])
  peers = [os.path.dirname(row["path"]) for row in rows
           if os.path.basename(row["path"]) == name and row["path"] != entry["path"]]
  tail = 2
  abbreviated = _recent_short_directory(directory, tail)
  # Show more of the path when otherwise identical names would be ambiguous.
  while any(other != directory and _recent_short_directory(other, tail) == abbreviated for other in peers):
    tail += 1
    abbreviated = _recent_short_directory(directory, tail)
  name += " | " + abbreviated
  if entry["kind"] == "mtz":
    name += " (MTZ; %d map%s)" % (len(entry["maps"]), "s" if len(entry["maps"]) != 1 else "")
  if not os.path.isfile(entry["path"]):
    name += " [missing]"
  return name


def _refresh_recent_files_menus():
  state = _recent_files_state()
  for category, menu in state["menus"].items():
    _clear_gio_menu(menu)
    rows = state["data"][category]
    if not rows:
      menu.append("No recent files yet", None)
    for index, entry in enumerate(rows):
      # Reuse action names: rebuilding a live menu must not leak app actions.
      # GMenu labels use underscores for mnemonics; double them only in the
      # menu presentation, never in the stored path or plain button labels.
      add_simple_action_to_menu(menu, _recent_file_label(entry, rows).replace("_", "__"),
        f"coot_trimmings_recent_{category}_{index}",
        lambda _action, _parameter, saved=entry: _open_recent_file(saved))
    for key, label, action in (
        ("refresh", "Refresh list", lambda category=category: _refresh_recent_files()),
        ("remove", "Remove entry...", lambda category=category: _prompt_remove_recent_file(category)),
        ("clear", "Clear list...", lambda category=category: _prompt_clear_recent_files(category))):
      add_simple_action_to_menu(menu, label, f"coot_trimmings_recent_{category}_{key}",
                                lambda _action, _parameter, action=action: action())


def _refresh_recent_files():
  _poll_recent_files()
  _save_recent_files()
  _refresh_recent_files_menus()


def _remove_recent_file(category, path):
  state = _recent_files_state()
  state["data"][category] = [entry for entry in state["data"][category] if entry["path"] != path]
  state["dirty"] = True
  _save_recent_files()
  _refresh_recent_files_menus()


def _prompt_remove_recent_file(category):
  rows = _recent_files_state()["data"][category]
  if not rows:
    info_dialog("No recent " + category + " are saved.")
    return
  action_button_dialog("Remove recent " + category,
    [(_recent_file_label(entry, rows), lambda path=entry["path"]: _remove_recent_file(category, path))
     for entry in rows], close_on_click=True)


def _prompt_clear_recent_files(category):
  def clear():
    state = _recent_files_state()
    state["data"][category] = []
    state["dirty"] = True
    _save_recent_files()
    _refresh_recent_files_menus()
  generic_confirm_dialog("Clear recent " + category,
    "Remove these history entries? The actual files will not be changed.", "Cancel", lambda: None,
    "Clear list", clear)


def _open_recent_file(entry):
  path = entry["path"]
  if not os.path.isfile(path):
    info_dialog("This recent file is missing or unavailable:\n\n" + path + "\n\nUse Remove entry to forget it.")
    _refresh_recent_files_menus()
    return None
  opened = []
  try:
    existing = {imol for imol in range(coot.graphics_n_molecules())
                if coot.is_valid_model_molecule(imol) or coot.is_valid_map_molecule(imol)}
    def valid_new(imol, model=False):
      valid = coot.is_valid_model_molecule if model else coot.is_valid_map_molecule
      return type(imol) is int and imol >= 0 and imol not in existing and valid(imol)
    if entry["kind"] == "model":
      imol = coot.read_pdb(path)
      if not valid_new(imol, model=True):
        raise RuntimeError("Coot could not read this model.")
      opened.append(imol)
    elif entry["kind"] == "map":
      imol = coot.read_ccp4_map(path, int(entry["difference"]))
      if not valid_new(imol):
        raise RuntimeError("Coot could not read this map.")
      opened.append(imol)
    else:
      for row in entry["maps"]:
        imol = coot.make_and_draw_map(path, row["f"], row["phi"], row["weight"],
                                      int(row["use_weights"]), int(row["difference"]))
        if not valid_new(imol):
          raise RuntimeError("Could not load map columns %s / %s; the MTZ may have changed." % (row["f"], row["phi"]))
        opened.append(imol)
    _poll_recent_files()
    add_status_bar_text("Opened recent file: " + path)
    return opened
  except Exception as error:
    # Only close molecules created by this operation, never existing models/maps.
    for imol in opened:
      coot.close_molecule(imol)
    info_dialog("Could not reopen recent file:\n\n" + path + "\n\n" + str(error))
    return None


def _install_recent_files_tracking(models_menu, maps_menu):
  state = _recent_files_state()
  state["menus"] = {"models": models_menu, "maps": maps_menu}
  glib = _coot_gui_repository_module("GLib")
  if state["timer"] is not None and glib is not None:
    glib.source_remove(state["timer"])
    state["timer"] = None
  _poll_recent_files()
  _refresh_recent_files_menus()
  if glib is not None:
    state["timer"] = glib.timeout_add_seconds(2, _poll_recent_files)


@lru_cache(maxsize=1)
def _ccp4_monomer_library_roots():
  roots = []
  candidates = [
    os.environ.get("CLIBD_MON"),
    os.path.join(os.environ.get("CCP4", ""), "lib", "data", "monomers") if os.environ.get("CCP4") else "",
  ]
  # Prefer the running Coot's data directory over another installed suite.
  # This also follows Homebrew upgrades without querying brew at runtime.
  package_data = getattr(coot, "get_pkgdatadir_py", None)
  if callable(package_data):
    try:
      data_dir = os.fspath(package_data())
      if data_dir:
        candidates.append(os.path.join(data_dir, "lib", "data", "monomers"))
    except Exception:
      pass
  for prefix in (os.environ.get("HOMEBREW_PREFIX"), "/opt/homebrew", "/usr/local", "/home/linuxbrew/.linuxbrew"):
    if prefix:
      candidates.append(os.path.join(prefix, "opt", "coot", "share", "coot", "lib", "data", "monomers"))
  candidates.append(
    "/Applications/CCP4/ccp4-9/lib/data/monomers",
  )
  seen = set()
  for candidate in candidates:
    if candidate and os.path.isdir(candidate):
      resolved = os.path.realpath(candidate)
      if resolved not in seen:
        seen.add(resolved)
        roots.append(candidate)
  return roots


@lru_cache(maxsize=None)
def _find_ccp4_monomer_file(monomer_code):
  normalized_code = str(monomer_code).strip().upper()
  if not normalized_code:
    return None
  for root in _ccp4_monomer_library_roots():
    candidate = os.path.join(root, normalized_code[0].lower(), normalized_code + ".cif")
    if os.path.isfile(candidate):
      return candidate
  return None


def _is_pointer_type_common_monomer(monomer_code):
  return bool(_pointer_type_for_monomer_code(monomer_code))


def _single_atom_pointer_type_from_cif(monomer_code):
  normalized_code = str(monomer_code).strip().upper()
  monomer_file = _find_ccp4_monomer_file(monomer_code)
  if not monomer_file:
    return None

  try:
    with open(monomer_file, "r", encoding="utf-8") as handle:
      lines = handle.readlines()
  except Exception:
    return None

  in_atom_loop = False
  atom_headers = []
  atom_rows = []

  for raw_line in lines:
    line = raw_line.strip()
    if not line or line.startswith("#"):
      continue
    if line == "loop_":
      in_atom_loop = False
      atom_headers = []
      continue
    if line.startswith("_chem_comp_atom."):
      in_atom_loop = True
      atom_headers.append(line)
      continue
    if in_atom_loop and atom_headers and line.startswith("_"):
      in_atom_loop = False
      atom_headers = []
      continue
    if in_atom_loop and atom_headers:
      atom_rows.append(line.split())
      if len(atom_rows) > 1:
        return None

  if len(atom_rows) != 1 or not atom_headers:
    return None

  try:
    type_symbol_index = atom_headers.index("_chem_comp_atom.type_symbol")
  except ValueError:
    return None

  atom_row = atom_rows[0]
  if type_symbol_index >= len(atom_row):
    return None

  type_symbol = atom_row[type_symbol_index].strip().upper()
  if not type_symbol:
    return None
  return normalized_code


def _pointer_type_for_monomer_code(monomer_code):
  normalized_code = str(monomer_code).strip().upper()
  if not normalized_code:
    return None
  explicit_pointer_type = COMMON_MONOMER_POINTER_TYPES.get(normalized_code)
  if explicit_pointer_type:
    return explicit_pointer_type
  return _single_atom_pointer_type_from_cif(normalized_code)


def _common_monomer_code_is_supported(monomer_code):
  normalized_code = str(monomer_code).strip().upper()
  if not normalized_code:
    return False
  if _pointer_type_for_monomer_code(normalized_code):
    return True
  return bool(_find_ccp4_monomer_file(normalized_code))


def _sanitize_common_monomer_favorite_name(name):
  allowed_punctuation = set(" -+/'.(),[]&:")
  cleaned = "".join(
    character
    for character in str(name)
    if character.isalnum() or character in allowed_punctuation
  )
  return " ".join(cleaned.split()).strip()


def _sanitize_common_monomer_code(code):
  return "".join(character for character in str(code).upper() if character.isalnum())


def _normalize_common_monomer_favorite_code(name, code):
  favorite_name = _sanitize_common_monomer_favorite_name(name)
  monomer_code = _sanitize_common_monomer_code(code)
  if not favorite_name or not monomer_code:
    return None
  return {"kind": "code", "name": favorite_name, "code": monomer_code}


def _normalize_common_monomer_cif_path(cif_path):
  normalized_path = os.path.abspath(os.path.expanduser(str(cif_path).strip()))
  if not normalized_path or not os.path.isfile(normalized_path):
    return None
  return normalized_path


@lru_cache(maxsize=None)
def _extract_cif_component_id(cif_path):
  try:
    with open(cif_path, "r", encoding="utf-8") as handle:
      for raw_line in handle:
        line = raw_line.strip()
        if not line or line.startswith("#"):
          continue
        if line.startswith("_chem_comp.id"):
          parts = line.split(maxsplit=1)
          if len(parts) == 2:
            comp_id = _sanitize_common_monomer_code(parts[1].strip().strip("'\""))
            if comp_id:
              return comp_id
        if line.startswith("data_comp_") and line != "data_comp_list":
          comp_id = _sanitize_common_monomer_code(line[len("data_comp_"):])
          if comp_id:
            return comp_id
  except Exception:
    return None
  return None


def _pointer_type_for_cif_path(cif_path, comp_id=None):
  normalized_path = _normalize_common_monomer_cif_path(cif_path)
  normalized_code = _sanitize_common_monomer_code(comp_id or "")
  if not normalized_path:
    return None

  try:
    with open(normalized_path, "r", encoding="utf-8") as handle:
      lines = handle.readlines()
  except Exception:
    return None

  in_atom_loop = False
  atom_headers = []
  atom_rows = []

  for raw_line in lines:
    line = raw_line.strip()
    if not line or line.startswith("#"):
      continue
    if line == "loop_":
      in_atom_loop = False
      atom_headers = []
      continue
    if line.startswith("_chem_comp_atom."):
      in_atom_loop = True
      atom_headers.append(line)
      continue
    if in_atom_loop and atom_headers and line.startswith("_"):
      in_atom_loop = False
      atom_headers = []
      continue
    if in_atom_loop and atom_headers:
      atom_rows.append(line.split())
      if len(atom_rows) > 1:
        return None

  if len(atom_rows) != 1 or not atom_headers:
    return None

  try:
    type_symbol_index = atom_headers.index("_chem_comp_atom.type_symbol")
  except ValueError:
    return None

  atom_row = atom_rows[0]
  if type_symbol_index >= len(atom_row):
    return None

  type_symbol = atom_row[type_symbol_index].strip().upper()
  if not type_symbol:
    return None
  return normalized_code or _extract_cif_component_id(normalized_path)


def _normalize_common_monomer_favorite_cif(name, cif_path, comp_id=None):
  favorite_name = _sanitize_common_monomer_favorite_name(name)
  normalized_path = _normalize_common_monomer_cif_path(cif_path)
  if not favorite_name or not normalized_path:
    return None
  normalized_comp_id = _sanitize_common_monomer_code(comp_id or "")
  if not normalized_comp_id:
    normalized_comp_id = _extract_cif_component_id(normalized_path)
  if not normalized_comp_id:
    return None
  return {
    "kind": "cif",
    "name": favorite_name,
    "cif_path": normalized_path,
    "comp_id": normalized_comp_id,
  }


def _validate_custom_cif_dictionary(cif_path, comp_id):
  normalized_path = _normalize_common_monomer_cif_path(cif_path)
  normalized_code = _sanitize_common_monomer_code(comp_id)
  if not normalized_path or not normalized_code:
    return None

  try:
    read_cif_dictionary(normalized_path)
  except Exception:
    traceback.print_exc()
    return None

  temp_imol = coot.get_monomer_from_dictionary(normalized_code, 1)
  if not valid_model_molecule_qm(temp_imol):
    return None

  try:
    return {
      "cif_path": normalized_path,
      "comp_id": normalized_code,
      "pointer_type": _pointer_type_for_cif_path(normalized_path, normalized_code),
    }
  finally:
    if valid_model_molecule_qm(temp_imol):
      coot.close_molecule(temp_imol)


def _managed_common_monomer_cif_path(source_cif_path, comp_id):
  normalized_path = _normalize_common_monomer_cif_path(source_cif_path)
  normalized_code = _sanitize_common_monomer_code(comp_id)
  if not normalized_path or not normalized_code:
    return None
  with open(normalized_path, "rb") as handle:
    digest = hashlib.sha1(handle.read()).hexdigest()[:12]
  return os.path.join(
    _coot_trimmings_dir(),
    f"{COMMON_MONOMER_FAVORITE_CIF_PREFIX}{normalized_code}_{digest}.cif",
  )


def _materialize_common_monomer_favorite_cif(name, cif_path):
  normalized = _normalize_common_monomer_favorite_cif(name, cif_path)
  if not normalized:
    return None
  validation = _validate_custom_cif_dictionary(normalized["cif_path"], normalized["comp_id"])
  if not validation:
    return None
  managed_cif_path = _managed_common_monomer_cif_path(normalized["cif_path"], normalized["comp_id"])
  if not managed_cif_path:
    return None
  os.makedirs(os.path.dirname(managed_cif_path), exist_ok=True)
  if os.path.abspath(normalized["cif_path"]) != os.path.abspath(managed_cif_path):
    shutil.copyfile(normalized["cif_path"], managed_cif_path)
  normalized["cif_path"] = managed_cif_path
  return normalized


def _favorite_entry_token(entry):
  if entry.get("kind") == "cif":
    return entry["comp_id"]
  return entry["code"]


def _favorite_menu_label(entry):
  if entry.get("kind") == "cif":
    return str(entry["name"]).strip()
  return f"{str(entry['name']).strip()} ({_favorite_entry_token(entry)})"


def _favorite_identity(entry):
  if entry.get("kind") == "cif":
    return ("cif", entry["cif_path"], entry["comp_id"])
  return ("code", entry["code"])


def _load_common_monomer_favorites():
  favorites_path = _common_monomer_favorites_path()
  if not os.path.isfile(favorites_path):
    return []

  try:
    with open(favorites_path, "r", encoding="utf-8") as handle:
      raw_favorites = json.load(handle)
  except Exception:
    print("coot_trimmings favorite-load failure:", favorites_path)
    traceback.print_exc()
    return []

  favorites = []
  seen_favorites = set()
  if isinstance(raw_favorites, list):
    for entry in raw_favorites:
      if not isinstance(entry, dict):
        continue
      if entry.get("kind") == "cif" or "cif_path" in entry:
        normalized = _normalize_common_monomer_favorite_cif(
          entry.get("name", ""),
          entry.get("cif_path", ""),
          entry.get("comp_id", ""),
        )
      else:
        normalized = _normalize_common_monomer_favorite_code(
          entry.get("name", ""),
          entry.get("code", ""),
        )
      if not normalized:
        continue
      favorite_identity = _favorite_identity(normalized)
      if favorite_identity in seen_favorites:
        continue
      if normalized["kind"] == "code":
        monomer_code = normalized["code"]
        if monomer_code not in COMMON_MONOMER_POINTER_TYPES and not _find_ccp4_monomer_file(monomer_code):
          continue
      if normalized["kind"] == "cif" and not os.path.isfile(normalized["cif_path"]):
        continue
      seen_favorites.add(favorite_identity)
      favorites.append(normalized)

  favorites.sort(key=lambda entry: _favorite_menu_label(entry).lower())
  return favorites


def _save_common_monomer_favorites(favorites):
  favorites_path = _common_monomer_favorites_path()
  if not favorites:
    if os.path.isfile(favorites_path):
      os.remove(favorites_path)
    return favorites_path
  os.makedirs(os.path.dirname(favorites_path), exist_ok=True)
  with open(favorites_path, "w", encoding="utf-8") as handle:
    json.dump(favorites, handle, indent=2, sort_keys=True)
  return favorites_path


def _clear_gio_menu(menu):
  if menu is None:
    return None
  if hasattr(menu, "remove_all"):
    menu.remove_all()
    return None
  if hasattr(menu, "get_n_items") and hasattr(menu, "remove"):
    while menu.get_n_items() > 0:
      menu.remove(0)
  return None


def refresh_common_monomer_favorites_menu():
  global COMMON_MONOMER_FAVORITES_MENU
  menu = COMMON_MONOMER_FAVORITES_MENU
  if menu is None:
    return None

  _clear_gio_menu(menu)
  add_simple_coot_menu_menuitem(menu, "Add favorite from monomer code...", lambda func: prompt_add_common_monomer_code_favorite())
  add_simple_coot_menu_menuitem(menu, "Add favorite from CIF...", lambda func: prompt_add_common_monomer_cif_favorite())
  add_simple_coot_menu_menuitem(menu, "Remove favorite...", lambda func: prompt_remove_common_monomer_favorite())

  for favorite in _load_common_monomer_favorites():
    add_simple_coot_menu_menuitem(
      menu,
      _favorite_menu_label(favorite),
      lambda func, favorite_entry=favorite: place_common_monomer_favorite(favorite_entry),
    )
  return None


def _upsert_common_monomer_code_favorite(name, code):
  normalized = _normalize_common_monomer_favorite_code(name, code)
  if not normalized:
    info_dialog(
      "Favorites require both a display name and a monomer code.\n\n"
      "Invalid characters are stripped automatically; if that leaves either field empty, the favorite is rejected."
    )
    return 0

  monomer_code = normalized["code"]
  if not _common_monomer_code_is_supported(monomer_code):
    info_dialog(
      f"{monomer_code} was not found in the CCP4 monomer library.\n\n"
      "Favorites are validated case-insensitively against the installed CCP4 monomer files,\n"
      "with typed-atom special cases handled for metal/ion monomers placed at the pointer."
    )
    return 0

  favorites = _load_common_monomer_favorites()
  favorite_identity = _favorite_identity(normalized)
  if favorite_identity in [_favorite_identity(entry) for entry in favorites]:
    favorites = [entry for entry in favorites if _favorite_identity(entry) != favorite_identity]
    status_text = f"Updated favorite: {_favorite_menu_label(normalized)}"
  else:
    status_text = f"Added favorite: {_favorite_menu_label(normalized)}"

  favorites.append(normalized)
  favorites.sort(key=lambda entry: _favorite_menu_label(entry).lower())
  try:
    favorites_path = _save_common_monomer_favorites(favorites)
  except Exception:
    info_dialog("Failed to save Common monomer favorites.")
    traceback.print_exc()
    return 0
  refresh_common_monomer_favorites_menu()
  add_status_bar_text(f"{status_text} [{favorites_path}]")
  return 1


def _upsert_common_monomer_cif_favorite(name, cif_path):
  normalized = _materialize_common_monomer_favorite_cif(name, cif_path)
  if not normalized:
    info_dialog(
      "CIF favorites require both a display name and a valid CIF dictionary.\n\n"
      "The file must exist, contain a readable component ID, and successfully build the monomer in Coot."
    )
    return 0

  favorites = _load_common_monomer_favorites()
  favorite_identity = _favorite_identity(normalized)
  if favorite_identity in [_favorite_identity(entry) for entry in favorites]:
    favorites = [entry for entry in favorites if _favorite_identity(entry) != favorite_identity]
    status_text = f"Updated favorite: {_favorite_menu_label(normalized)}"
  else:
    status_text = f"Added favorite: {_favorite_menu_label(normalized)}"

  favorites.append(normalized)
  favorites.sort(key=lambda entry: _favorite_menu_label(entry).lower())
  try:
    favorites_path = _save_common_monomer_favorites(favorites)
  except Exception:
    info_dialog("Failed to save Common monomer favorites.")
    traceback.print_exc()
    return 0
  refresh_common_monomer_favorites_menu()
  add_status_bar_text(f"{status_text} [{favorites_path}]")
  return 1


def prompt_add_common_monomer_code_favorite():
  generic_double_entry(
    "Add Common Monomer Favorite",
    "Favorite name",
    "",
    "Monomer code",
    "",
    "Add favorite",
    _upsert_common_monomer_code_favorite,
  )
  return 1


def prompt_add_common_monomer_cif_favorite():
  generic_double_entry_with_file_browse(
    "Add Common Monomer Favorite",
    "Favorite name",
    "",
    "CIF file path",
    "",
    "Choose CIF dictionary",
    "Add favorite",
    _upsert_common_monomer_cif_favorite,
  )
  return 1


def _remove_common_monomer_favorite(favorite_entry):
  favorites = _load_common_monomer_favorites()
  favorite_identity = _favorite_identity(favorite_entry)
  new_favorites = [entry for entry in favorites if _favorite_identity(entry) != favorite_identity]
  if len(new_favorites) == len(favorites):
    info_dialog(f"{_favorite_menu_label(favorite_entry)} is not currently in Common monomer favorites.")
    return 0
  try:
    _save_common_monomer_favorites(new_favorites)
  except Exception:
    info_dialog("Failed to save Common monomer favorites.")
    traceback.print_exc()
    return 0
  if favorite_entry.get("kind") == "cif":
    managed_cif_path = favorite_entry.get("cif_path", "")
    if (
      isinstance(managed_cif_path, str)
      and managed_cif_path
      and os.path.dirname(os.path.abspath(managed_cif_path)) == _coot_trimmings_dir()
      and os.path.basename(managed_cif_path).startswith(COMMON_MONOMER_FAVORITE_CIF_PREFIX)
      and not any(
        entry.get("kind") == "cif" and entry.get("cif_path") == managed_cif_path
        for entry in new_favorites
      )
    ):
      try:
        if os.path.isfile(managed_cif_path):
          os.remove(managed_cif_path)
      except Exception:
        traceback.print_exc()
  refresh_common_monomer_favorites_menu()
  add_status_bar_text(f"Removed favorite: {_favorite_menu_label(favorite_entry)}")
  return 1


def prompt_remove_common_monomer_favorite():
  favorites = _load_common_monomer_favorites()
  if not favorites:
    info_dialog("No Common monomer favorites are currently saved.")
    return 0

  action_button_dialog(
    "Remove Common Monomer Favorite",
    [
      (
        _favorite_menu_label(entry),
        lambda favorite_entry=entry: _remove_common_monomer_favorite(favorite_entry),
      )
      for entry in favorites
    ],
    close_on_click=True,
  )
  return 1


def place_common_monomer_favorite(favorite_entry):
  if favorite_entry.get("kind") == "cif":
    return place_custom_cif_monomer(favorite_entry["cif_path"], favorite_entry["comp_id"])
  return place_common_monomer(favorite_entry["code"])


def _common_monomer_target_molecule():
  residue = active_residue()
  if residue:
    return residue[0]
  try:
    mol_id = go_to_atom_molecule_number()
    if valid_model_molecule_qm(mol_id):
      return mol_id
  except Exception:
    pass
  for mol_id in model_molecule_list():
    if valid_model_molecule_qm(mol_id):
      return mol_id
  return None


def place_custom_cif_monomer(cif_path, comp_id=None):
  normalized_path = _normalize_common_monomer_cif_path(cif_path)
  normalized_code = _sanitize_common_monomer_code(comp_id or _extract_cif_component_id(cif_path) or "")
  if not normalized_path or not normalized_code:
    info_dialog("This favorite references an invalid CIF file or component ID.")
    return 0

  pointer_type = _pointer_type_for_cif_path(normalized_path, normalized_code)
  if pointer_type:
    target_mol_id = _common_monomer_target_molecule()
    if target_mol_id is None:
      info_dialog(
        f"{normalized_code} could not be added because there is no target model molecule.\n\n"
        "Single-atom monomers such as metals and halides are inserted into an existing model.\n"
        "Load a model first, or make sure a model molecule is active, then try again."
      )
      return 0
    try:
      set_pointer_atom_molecule(target_mol_id)
    except Exception:
      pass
    try:
      set_go_to_atom_molecule(target_mol_id)
    except Exception:
      pass
    place_typed_atom_at_pointer(pointer_type)
    return target_mol_id

  read_cif_dictionary(normalized_path)
  monomer_imol = coot.get_monomer_from_dictionary(normalized_code, 1)
  if not valid_model_molecule_qm(monomer_imol):
    monomer_imol = coot.get_monomer(normalized_code)
  if not valid_model_molecule_qm(monomer_imol):
    info_dialog(f"Failed to build {normalized_code} from {normalized_path}.")
    return 0
  delete_hydrogens(monomer_imol)
  return monomer_imol


def place_common_monomer(monomer_code):
  pointer_type = _pointer_type_for_monomer_code(monomer_code)
  if pointer_type:
    target_mol_id = _common_monomer_target_molecule()
    if target_mol_id is None:
      info_dialog(
        f"{str(monomer_code).strip().upper()} could not be added because there is no target model molecule.\n\n"
        "Single-atom monomers such as metals and halides are inserted into an existing model.\n"
        "Load a model first, or make sure a model molecule is active, then try again."
      )
      return 0
    try:
      set_pointer_atom_molecule(target_mol_id)
    except Exception:
      pass
    try:
      set_go_to_atom_molecule(target_mol_id)
    except Exception:
      pass
    place_typed_atom_at_pointer(pointer_type)
    return target_mol_id
  return get_monomer_no_H(monomer_code)


def add_common_monomer_menu_entries(menu, entries):
  if menu is None:
    return None
  gio_module = _gio_module()
  if gio_module is None:
    return None
  for label, value in entries:
    if isinstance(value, list):
      submenu = gio_module.Menu.new()
      menu.append_submenu(label, submenu)
      add_common_monomer_menu_entries(submenu, value)
    else:
      add_simple_coot_menu_menuitem(
        menu,
        label,
        lambda func, monomer_code=value: place_common_monomer(monomer_code),
      )
  return None


COORDINATION_LINK_MENU = [
  ("O donors", [
    ("Mg-O coordination", {"metal_codes": ("MG",), "metal_label": "Mg", "donor_elements": ("O",), "donor_label": "oxygen donor", "distance": 2.10, "range_text": "expected range ~2.0-2.2 A"}),
    ("Na-O coordination", {"metal_codes": ("NA",), "metal_label": "Na", "donor_elements": ("O",), "donor_label": "oxygen donor", "distance": 2.41, "range_text": "expected range ~2.3-2.6 A"}),
    ("K-O coordination", {"metal_codes": ("K",), "metal_label": "K", "donor_elements": ("O",), "donor_label": "oxygen donor", "distance": 2.80, "range_text": "expected range ~2.7-3.2 A"}),
    ("Ca-O coordination", {"metal_codes": ("CA",), "metal_label": "Ca", "donor_elements": ("O",), "donor_label": "oxygen donor", "distance": 2.40, "range_text": "expected range ~2.3-2.6 A"}),
    ("Zn-O coordination", {"metal_codes": ("ZN",), "metal_label": "Zn", "donor_elements": ("O",), "donor_label": "oxygen donor", "distance": 2.00, "range_text": "expected range ~1.9-2.2 A"}),
    ("Fe-O coordination", {"metal_codes": ("FE",), "metal_label": "Fe", "donor_elements": ("O",), "donor_label": "oxygen donor", "distance": 2.05, "range_text": "expected range ~2.0-2.2 A"}),
  ]),
  ("N donors", [
    ("Zn-N coordination", {"metal_codes": ("ZN",), "metal_label": "Zn", "donor_elements": ("N",), "donor_label": "nitrogen donor", "distance": 2.058, "range_text": "expected range ~1.9-2.2 A"}),
    ("Fe-N coordination", {"metal_codes": ("FE",), "metal_label": "Fe", "donor_elements": ("N",), "donor_label": "nitrogen donor", "distance": 2.058, "range_text": "expected range ~2.0-2.24 A"}),
  ]),
  ("S donors", [
    ("Zn-S coordination", {"metal_codes": ("ZN",), "metal_label": "Zn", "donor_elements": ("S",), "donor_label": "sulfur donor", "distance": 2.340, "range_text": "expected range ~2.2-2.4 A"}),
    ("Fe-S coordination", {"metal_codes": ("FE",), "metal_label": "Fe", "donor_elements": ("S",), "donor_label": "sulfur donor", "distance": 2.260, "range_text": "expected range ~2.27-2.36 A"}),
  ]),
]


# Reference fallbacks, NOT oxidation-state/geometry-specific PDB ideal distances.
# Distances: https://github.com/pemsley/coot/tree/main/data/metal
# Counts and donor-element restrictions: cctbx/mmtbx/ions/ion_parameters.cif
# https://github.com/cctbx/cctbx_project/blob/master/mmtbx/ions/ion_parameters.cif
# That file cites Dokmanic et al. (2008), Zheng et al. (2008), Harding (2001),
# and Rulisek & Vondrasek (1998). Its bond-valence parameters are NOT lengths.
# Embedded so neither Phenix nor a network connection is needed at runtime.
METAL_LINK_REFERENCE_PROFILES = {
  "MN": ((4, 6), {"O": 1.9159, "N": 1.9329}, "Mn(II)"),
  "CO": ((5, 6), {"O": 1.9221, "N": 1.8284}, "Co(II)"),
  "NI": ((4, 6), {"O": 2.0544, "N": 2.1043, "S": 2.1638}, "Ni(II)"),
  "CU": ((3, 4), {"O": 1.9888, "N": 2.5317, "S": 2.1835}, "Cu(II)"),
  "CD": ((3, 7), {"O": 2.3494, "N": 2.3164, "S": 2.9976}, "Cd(II)"),
}


# BEGIN GENERATED METALPDB PROFILES
# Generated by tools/build_metalpdb_profiles.py; provenance: tools/metalpdb_reference/profiles.json.
# Observed medians, NOT ideal restraints. No resolution/oxidation-state filter.
# MetalPDB truncates contacts at 3 A; p90 >=2.9 A excluded. Counts are not chemical limits.
# Per donor class: most-supported geometry/CN stratum, not inferred geometry of this site.
# Values: distance, geometry, CN, PDB count, p10, p90 (per-PDB median distances).
METALPDB_LINK_PROFILES = {
  'LI': ((4, 4), {
    'carboxylate-O': (1.919, 'tetrahedron (regular)', 4, 5, 1.822, 1.988),
    'water-O': (1.957, 'tetrahedron (regular)', 4, 5, 1.877, 1.973),
  }),
  'SR': ((4, 8), {
    'carboxylate-O': (2.64, 'trigonal prism, square-face bicapped (distorted)', 8, 5, 2.6, 2.731),
    'water-O': (2.742, 'irregular (n/a)', 4, 7, 2.634, 2.763),
  }),
  'HG': ((2, 4), {
    'backbone-O': (2.755, 'irregular (n/a)', 4, 6, 2.663, 2.896),
    'carboxylate-O': (2.637, 'irregular (n/a)', 3, 6, 2.443, 2.77),
    'cysteine-S': (2.397, 'linear (regular)', 2, 45, 2.298, 2.614),
  }),
  'PR': ((3, 5), {
    'carboxylate-O': (2.729, 'irregular (n/a)', 3, 8, 2.545, 2.828),
    'water-O': (2.615, 'irregular (n/a)', 3, 5, 2.504, 2.892),
  }),
  'SM': ((4, 7), {
    'carboxylate-O': (2.503, 'irregular (n/a)', 4, 7, 2.337, 2.675),
    'water-O': (2.524, 'square antiprism with a vacancy (distorted)', 7, 5, 2.23, 2.664),
  }),
  'GD': ((2, 9), {
    'carboxylate-O': (2.641, 'irregular (n/a)', 2, 5, 2.503, 2.727),
    'water-O': (2.452, 'irregular (n/a)', 9, 6, 2.322, 2.625),
  }),
  'YB': ((3, 5), {
    'carboxylate-O': (2.405, 'irregular (n/a)', 5, 7, 2.304, 2.612),
    'water-O': (2.506, 'irregular (n/a)', 3, 5, 2.22, 2.65),
  }),
  'MN': ((2, 8), {
    'amide-O': (2.177, 'octahedron (regular)', 6, 24, 2.046, 2.348),
    'backbone-O': (2.195, 'octahedron (regular)', 6, 49, 2.101, 2.317),
    'carboxylate-O': (2.171, 'octahedron (regular)', 6, 173, 2.082, 2.305),
    'histidine-N': (2.254, 'octahedron (regular)', 6, 71, 2.18, 2.371),
    'hydroxyl-O': (2.227, 'octahedron (regular)', 6, 20, 2.131, 2.372),
    'water-O': (2.233, 'octahedron (regular)', 6, 200, 2.109, 2.367),
  }),
  'CO': ((2, 6), {
    'backbone-O': (2.122, 'octahedron (regular)', 6, 8, 2.026, 2.206),
    'carboxylate-O': (2.11, 'octahedron (regular)', 6, 41, 1.978, 2.222),
    'cysteine-S': (2.294, 'tetrahedron (regular)', 4, 5, 2.138, 2.44),
    'histidine-N': (2.171, 'octahedron (regular)', 6, 43, 2.047, 2.266),
    'hydroxyl-O': (2.19, 'octahedron (regular)', 6, 7, 2.005, 2.421),
    'water-O': (2.187, 'octahedron (regular)', 6, 51, 2.053, 2.343),
  }),
  'NI': ((2, 7), {
    'amide-O': (2.291, 'octahedron (regular)', 6, 8, 2.053, 2.481),
    'backbone-O': (2.12, 'octahedron (regular)', 6, 16, 2.051, 2.295),
    'carboxylate-O': (2.125, 'octahedron (regular)', 6, 34, 2.027, 2.239),
    'cysteine-S': (2.23, 'square plane (regular)', 4, 7, 2.077, 2.338),
    'histidine-N': (2.133, 'octahedron (regular)', 6, 74, 2.046, 2.258),
    'water-O': (2.176, 'octahedron (regular)', 6, 67, 2.031, 2.388),
  }),
  'CU': ((2, 6), {
    'carboxylate-O': (2.401, 'irregular (n/a)', 3, 14, 2.113, 2.513),
    'cysteine-S': (2.222, 'trigonal plane (regular)', 3, 13, 2.186, 2.27),
    'histidine-N': (2.037, 'irregular (n/a)', 5, 32, 1.941, 2.21),
    'methionine-S': (2.304, 'tetrahedron (regular)', 4, 7, 2.214, 2.43),
    'water-O': (1.838, 'linear (regular)', 2, 6, 1.797, 2.041),
  }),
  'CD': ((2, 7), {
    'amide-O': (2.362, 'octahedron (distorted)', 6, 6, 2.244, 2.482),
    'backbone-O': (2.626, 'irregular (n/a)', 3, 17, 2.383, 2.828),
    'carboxylate-O': (2.466, 'irregular (n/a)', 3, 150, 2.302, 2.726),
    'cysteine-S': (2.544, 'tetrahedron (regular)', 4, 22, 2.455, 2.6),
    'histidine-N': (2.326, 'irregular (n/a)', 5, 54, 2.191, 2.554),
    'hydroxyl-O': (2.443, 'octahedron (regular)', 6, 5, 2.37, 2.565),
    'water-O': (2.378, 'octahedron (regular)', 6, 40, 2.25, 2.514),
  }),
}
# END GENERATED METALPDB PROFILES

def _metal_donor_class(residue, atom, element):
  """Match the donor classes used by the offline MetalPDB profile builder."""
  if element == "O":
    if residue in {"HOH", "WAT", "DOD"}:
      return "water-O"
    if (residue == "ASP" and atom in {"OD1", "OD2"}) or (residue == "GLU" and atom in {"OE1", "OE2"}):
      return "carboxylate-O"
    if (residue == "ASN" and atom == "OD1") or (residue == "GLN" and atom == "OE1"):
      return "amide-O"
    if (residue, atom) in {("SER", "OG"), ("THR", "OG1"), ("TYR", "OH")}:
      return "hydroxyl-O"
    if residue in set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split()) and atom == "O":
      return "backbone-O"
  if residue == "HIS" and atom in {"ND1", "NE2"} and element == "N":
    return "histidine-N"
  if residue == "CYS" and atom == "SG" and element == "S":
    return "cysteine-S"
  if residue == "MET" and atom == "SD" and element == "S":
    return "methionine-S"
  return None


# Shannon (1976), Acta Cryst. A32, 751-767, doi:10.1107/S0567739476001551.
# Ionic radii verified against pymatgen v2024.5.1/core/periodic_table.json.
# CN VIII except Ag(I), CN IV. Ligand offsets below are deliberately rough
# size-based heuristics, not tabulated neutral-ligand bond lengths.
METAL_LINK_ESTIMATE_RADII = {
  "RB": (1.61, "Rb(I), radius CN VIII"), "CS": (1.74, "Cs(I), radius CN VIII"),
  "BA": (1.42, "Ba(II), radius CN VIII"), "AG": (1.00, "Ag(I), radius CN IV"),
  "LA": (1.160, "La(III), radius CN VIII"), "CE": (1.143, "Ce(III), radius CN VIII"),
  "ND": (1.109, "Nd(III), radius CN VIII"), "PM": (1.093, "Pm(III), radius CN VIII"),
  "EU": (1.066, "Eu(III), radius CN VIII"), "TB": (1.040, "Tb(III), radius CN VIII"),
  "DY": (1.027, "Dy(III), radius CN VIII"), "HO": (1.015, "Ho(III), radius CN VIII"),
  "ER": (1.004, "Er(III), radius CN VIII"), "TM": (0.994, "Tm(III), radius CN VIII"),
  "LU": (0.977, "Lu(III), radius CN VIII"),
}


def _metal_link_candidates(imol, metal_spec, radius, estimate=None, candidate_specs=None):
  """Return chemistry-filtered candidates; observed distances only rank the search."""
  profiles = {"MG": (6, 6), "NA": (4, 8), "K": (6, 9),
              "CA": (6, 8), "ZN": (4, 6), "FE": (4, 6)}
  profiles.update({metal: values[0] for metal, values in METAL_LINK_REFERENCE_PROFILES.items()})
  for metal, values in METALPDB_LINK_PROFILES.items():
    profiles.setdefault(metal, values[0])
  atoms = residue_info_py(imol, *metal_spec[:3]) or []
  ion = [a for a in atoms if list(a[0][:2]) == metal_spec[3:]]
  if len(ion) != 1:
    raise ValueError("Select a single metal atom in the active residue.")
  element = str(ion[0][1][2]).strip().upper()
  if estimate is not None:
    if element not in METAL_LINK_ESTIMATE_RADII or estimate.get("element") != element:
      raise ValueError("The estimate does not match the active metal.")
    if any(not math.isfinite(estimate[key]) or not 1 <= estimate[key] <= 4.5 for key in ("O", "N")):
      raise ValueError("Estimated distances must be between 1 and 4.5 A.")
    profiles[element] = (0, 0)  # No inferred coordination-count limits.
  if element not in profiles:
    raise ValueError("Supported metals: " + ", ".join(sorted(profiles)) + ".")
  targets = {}
  for _, entries in COORDINATION_LINK_MENU:
    for _, entry in entries:
      if element in entry["metal_codes"]:
        targets.update({donor: entry["distance"] for donor in entry["donor_elements"]})
  if element in METAL_LINK_REFERENCE_PROFILES:
    targets = METAL_LINK_REFERENCE_PROFILES[element][1]
  empirical = METALPDB_LINK_PROFILES.get(element)
  if empirical:
    targets = {role.rsplit("-", 1)[1]: None for role in empirical[1]}
  if estimate is not None:
    targets = {key: estimate[key] for key in ("O", "N")}
    empirical = None
  if not callable(getattr(coot, "link_info_py", None)):
    raise ValueError("This Coot build cannot inspect existing links safely.")
  existing = set()
  for link in coot.link_info_py(imol) or []:
    a, b = tuple(link[1][-5:]), tuple(link[2][-5:])
    if a == tuple(metal_spec):
      existing.add(b)
    elif b == tuple(metal_spec):
      existing.add(a)
  protein = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split())
  candidates = []
  if candidate_specs is None:
    candidate_specs = _residue_specs_near_point(imol, ion[0][2], radius)
  for spec in candidate_specs:
    resname = residue_name(imol, *spec).strip().upper()
    ligand = resname not in protein and resname not in ("HOH", "WAT", "DOD")
    for atom in residue_info_py(imol, *spec) or []:
      name, alt = atom[0][:2]
      donor = str(atom[1][2]).strip().upper()
      if donor not in targets or atom[1][0] <= 0:
        continue
      # Do not mix alternate conformers, or choose one silently for a blank ion.
      if alt and alt != metal_spec[4]:
        continue
      atom_name = name.strip()
      reference = None
      proxy_role = None
      if empirical:
        role = _metal_donor_class(resname, atom_name, donor)
        reference = empirical[1].get(role)
        if reference is None and ligand:
          # A donor-element proxy permits manual ligand review without claiming
          # that protein/water statistics describe this ligand's chemistry.
          alternatives = [(key, value) for key, value in empirical[1].items()
                          if key.rsplit("-", 1)[1] == donor]
          if alternatives:
            proxy_role, reference = max(alternatives, key=lambda item: item[1][3])
        if reference is None:
          continue
      unlikely_n = False
      if donor == "N":
        plausible_n = ligand or (resname == "HIS" and atom_name in ("ND1", "NE2"))
        if estimate is not None:
          plausible_n = plausible_n or resname not in protein
          unlikely_n = not plausible_n
          if unlikely_n and not estimate.get("unlikely_n", False):
            continue
        elif not plausible_n:
          continue
      if donor == "S" and not (ligand or (resname == "CYS" and atom_name == "SG") or
                               (resname == "MET" and atom_name == "SD")):
        continue
      observed = math.sqrt(sum((float(x) - float(y)) ** 2
                               for x, y in zip(atom[2], ion[0][2])))
      donor_spec = list(spec) + [name, alt]
      if not math.isfinite(observed) or observed > radius or tuple(donor_spec) in existing:
        continue
      chemistry = "oxygen donor" if donor == "O" else "check protonation/oxidation state"
      if resname not in protein and resname not in ("HOH", "WAT", "DOD"):
        chemistry = "nonstandard residue: verify donor chemistry manually"
      target = targets[donor]
      if reference:
        target, geometry, cn, count, low, high = reference
        chemistry = (f"MetalPDB {proxy_role or role}; reference: {geometry}, CN {cn}. "
                     f"{count} PDB entries; central observed range {low:.3f}-{high:.3f} A. "
                     "Reference geometry is not inferred for the current site.")
        if count < 10:
          chemistry += " Limited evidence (fewer than 10 PDB entries)."
        if proxy_role:
          chemistry = "LOW CONFIDENCE: donor-element proxy, NOT ligand-specific data. " + chemistry
      if estimate is not None:
        chemistry = "LOW CONFIDENCE: editable size-based estimate; verify charge, donor chemistry and geometry."
        if unlikely_n:
          chemistry += " UNLIKELY N DONOR: amide/protonated nitrogen generally does not coordinate."
        elif donor == "N" and resname not in protein:
          chemistry += " Unclassified ligand N: bonding and protonation have not been checked."
      if ligand:
        chemistry += (" Ligand donor: verify bonding, protonation and chelation geometry; "
                      "nearby atoms are candidates, not an automatic chemical assignment.")
      candidates.append((donor_spec, target, observed, resname, chemistry))
  candidates.sort(key=lambda item: item[2])
  return element, profiles[element], len(existing), candidates


def _preview_metal_link(imol, metal_spec, donor_spec):
  """Own only this preview's drawing/text handles; leave atom labels untouched."""
  required = ("new_generic_object_number", "to_generic_object_add_dashed_line",
              "set_display_generic_object", "close_generic_object", "place_text", "remove_text")
  if any(not callable(getattr(coot, name, None)) for name in required):
    raise ValueError("This Coot build lacks the temporary link preview APIs.")
  positions = []
  for spec in (metal_spec, donor_spec):
    atoms = [a for a in residue_info_py(imol, *spec[:3]) or [] if list(a[0][:2]) == spec[3:]]
    if len(atoms) != 1:
      raise ValueError("A preview atom is no longer uniquely available.")
    positions.append(atoms[0][2])
  obj = coot.new_generic_object_number("Proposed metal coordination")
  labels = []
  closed = False
  def show_labels(enabled):
    for handle in reversed(labels):
      coot.remove_text(handle)
    labels.clear()
    if enabled and not closed:
      for spec, position in zip((metal_spec, donor_spec), positions):
        label = f"{spec[0]}:{spec[1]}{spec[2]} {spec[3].strip()}" + (f" alt {spec[4]}" if spec[4] else "")
        labels.append(coot.place_text(label, position[0] + 0.3, position[1] + 0.3, position[2], 14))
  def cleanup():
    nonlocal closed
    if closed:
      return
    closed = True
    for handle in reversed(labels):
      try:
        coot.remove_text(handle)
      except Exception as error:
        print(f"Could not remove preview label: {error}")
    coot.close_generic_object(obj)
  try:
    coot.to_generic_object_add_dashed_line(obj, "orange", 1, 5.0,
                                          *(list(positions[0]) + list(positions[1])))
    coot.set_display_generic_object(obj, 1)
    show_labels(True)
  except Exception:
    cleanup()
    raise
  return cleanup, show_labels


def _metal_link_review_panel(title, message, discard, keep, preview, state):
  """Non-modal review: the main Coot canvas remains interactive."""
  cleanup, show_labels = preview
  state.update(cleanup=cleanup, show_labels=show_labels, discard=discard, keep=keep)
  if "window" in state:
    state["window"].set_title(title)
    state["label"].set_text(message)
    show_labels(state["labels_button"].get_active())
    return
  window = Gtk.Window()
  window.set_title(title)
  window.set_modal(False)
  window.set_default_size(420, -1)
  box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=10)
  for side in ("start", "end", "top", "bottom"):
    getattr(box, "set_margin_" + side)(12)
  label = Gtk.Label(label=message)
  label.set_wrap(True)
  label.set_xalign(0)
  box.append(label)
  labels_button = Gtk.CheckButton(label="Show atom labels")
  labels_button.set_active(True)
  labels_button.connect("toggled", lambda button: state["show_labels"](button.get_active()))
  box.append(labels_button)
  buttons = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
  state.update(window=window, label=label, labels_button=labels_button, busy=False, closed=False)
  def finish(action=None):
    if state["closed"] or state["busy"]:
      return False
    state["busy"] = True
    try:
      state["cleanup"]()
      if action is None:
        state["closed"] = True
        window.destroy()
      else:
        state[action]()
    finally:
      state["busy"] = False
    return False
  for text, action in (("Finish", "finish_review"), ("Discard", "discard"), ("Keep link", "keep")):
    button = Gtk.Button(label=text)
    button.connect("clicked", lambda *_, action=action: finish(action))
    buttons.append(button)
  box.append(buttons)
  window.set_child(box)
  window.connect("close-request", lambda *_: finish())
  window.present()


def _metal_link_estimate_dialog(element, search):
  radius, assumption = METAL_LINK_ESTIMATE_RADII[element]
  window = Gtk.Window()
  window.set_title("Estimated metal coordination")
  box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=8)
  for side in ("start", "end", "top", "bottom"):
    getattr(box, "set_margin_" + side)(12)
  note = Gtk.Label(label=f"LOW CONFIDENCE: {assumption}.\n"
                   "Starting estimates: Shannon metal radius + 1.40 A (O) / 1.50 A (N).\n"
                   "Ligand offsets are heuristics, not measured bond lengths. Edit for your chemistry.\n"
                   "No coordination-number limits are inferred. All links require review.")
  note.set_wrap(True)
  box.append(note)
  entries = []
  for label, default in (("Search radius (A)", 5.0), ("Estimated metal-O (A)", radius + 1.4),
                         ("Estimated metal-N (A)", radius + 1.5)):
    box.append(Gtk.Label(label=label))
    entry = Gtk.Entry()
    entry.set_text(f"{default:.3f}")
    entries.append(entry)
    box.append(entry)
  unlikely = Gtk.CheckButton(label="Show unlikely N donors (amide/protonated nitrogens)")
  box.append(unlikely)
  def submit(*_args):
    try:
      estimate = {"element": element, "O": float(entries[1].get_text()),
                  "N": float(entries[2].get_text()), "unlikely_n": unlikely.get_active()}
      result = search(entries[0].get_text(), estimate)
    except ValueError:
      info_dialog("Please enter numeric distances.")
      return
    if result != 0:
      window.destroy()
  button = Gtk.Button(label="Review candidates")
  button.connect("clicked", submit)
  box.append(button)
  window.set_child(box)
  window.present()


def auto_metal_links_current_residue():
  active = _active_residue_or_status()
  if not active:
    return
  imol = active[0]
  metal_spec = list(active[1:6])

  def search(text, estimate=None):
    try:
      radius = float(text)
      if not math.isfinite(radius) or not 0 < radius <= 10:
        raise ValueError("Search radius must be greater than 0 and at most 10 A.")
      element, limits, existing_count, candidates = _metal_link_candidates(imol, metal_spec, radius, estimate)
    except (ValueError, TypeError, IndexError) as error:
      info_dialog(str(error))
      return 0
    if not candidates:
      info_dialog("No new compatible donors found. Existing links and incompatible alternate conformers are excluded.")
      return 0
    if estimate is not None:
      target_label = "editable estimate"
      source_note = f"LOW CONFIDENCE: {METAL_LINK_ESTIMATE_RADII[element][1]}. No expected coordination count inferred. "
    elif element in METALPDB_LINK_PROFILES:
      target_label = "MetalPDB observed median"
      source_note = ("Observed reference, not an ideal restraint. MetalPDB omits contacts beyond 3 A. No resolution or oxidation-state filter; "
                     "PDB entries may be homologous or restrained. Ligand proxies are explicitly flagged; "
                     "donor elements without reference data are omitted. ")
      if element in METAL_LINK_REFERENCE_PROFILES:
        source_note += f"Count guidance assumes {METAL_LINK_REFERENCE_PROFILES[element][2]}. "
      else:
        source_note += "Count range is observed in the retained reference strata, not a chemical limit. "
    elif element in METAL_LINK_REFERENCE_PROFILES:
      assumed_state = METAL_LINK_REFERENCE_PROFILES[element][2]
      target_label = "Coot reference distance"
      source_note = (f"Count guidance assumes {assumed_state}; oxidation state is not inferred. "
                     "Coot table distances are generic fallbacks, not validated PDB ideals; "
                     "verify the chemistry before keeping links. ")
    else:
      target_label = "chemistry target"
      source_note = ""
    selected = []
    empirical_counts_only = estimate is not None or (element in METALPDB_LINK_PROFILES and element not in METAL_LINK_REFERENCE_PROFILES)
    count_label = "Reference contact counts (not full coordination)" if empirical_counts_only else "Typical coordination"
    count_text = "Expected coordination: not inferred" if estimate is not None else f"{count_label}: {limits[0]}-{limits[1]}"
    panel_state = {}
    def close_panel():
      if "window" in panel_state and not panel_state.get("closed"):
        panel_state["closed"] = True
        panel_state["cleanup"]()
        panel_state["window"].destroy()

    def review(index):
      if index < len(candidates):
        spec, target, observed, resname, chemistry = candidates[index]
        try:
          preview = _preview_metal_link(imol, metal_spec, spec)
          cleanup = preview[0]
        except Exception as error:
          close_panel()
          info_dialog(f"Preview stopped; previously kept links remain.\n{error}")
          return
        def discard():
          cleanup()
          review(index + 1)
        def keep():
          cleanup()
          try:
            if not valid_model_molecule_qm(imol):
              raise ValueError("Target model is no longer available.")
            _, _, _, fresh = _metal_link_candidates(imol, metal_spec, radius, estimate,
                                                    candidate_specs=[tuple(spec[:3])])
            if tuple(spec) not in {tuple(item[0]) for item in fresh}:
              raise ValueError("The candidate changed or is already linked. Please restart the search.")
            coot.make_link_py(imol, metal_spec, spec, "dummy", target)
          except Exception as error:
            close_panel()
            info_dialog(f"Link creation stopped; previously kept links remain.\n{error}")
            return
          selected.append((spec, target))
          review(index + 1)
        warning = "\nWARNING: keeping this donor exceeds the typical coordination range." if not empirical_counts_only and existing_count + len(selected) + 1 > limits[1] else ""
        _metal_link_review_panel(
          f"{element} coordination candidate {index + 1}/{len(candidates)}",
          f"{element} {metal_spec[0]}:{metal_spec[1]}{metal_spec[2]} {metal_spec[3].strip()}\n"
          f"Donor: {resname} {spec[0]}:{spec[1]}{spec[2]} {spec[3].strip()}\n"
          f"Observed: {observed:.2f} A; {target_label}: {target:.3f} A\n{chemistry}\n{source_note}\n"
          f"Existing links: {existing_count}; kept: {len(selected)}. {count_text}.{warning}\n"
          "Dashed orange line is a temporary preview. You can rotate and zoom in Coot. Keep creates the real link immediately. "
          "Finish skips the remaining candidates and keeps accepted links, without accepting the current preview. "
          "Coot may use built-in refinement distances.",
          discard, keep, preview, panel_state)
        return
      close_panel()
      if not selected:
        info_dialog("No donors selected; no links added.")
        return
      total = existing_count + len(selected)
      warning = "\nWARNING: outside the typical coordination range." if not empirical_counts_only and not limits[0] <= total <= limits[1] else ""
      message = f"Submitted {len(selected)} {element} links; total coordination {total} ({count_text}).{warning}"
      print(message)
      add_status_bar_text(message)
      info_dialog(message)
    panel_state["finish_review"] = lambda: review(len(candidates))
    review(0)
    return None
  atoms = [a for a in residue_info_py(imol, *metal_spec[:3]) or [] if list(a[0][:2]) == metal_spec[3:]]
  element = str(atoms[0][1][2]).strip().upper() if len(atoms) == 1 else ""
  if element in METAL_LINK_ESTIMATE_RADII:
    _metal_link_estimate_dialog(element, search)
  else:
    generic_single_entry("Search radius around active metal (A)", "5.0", "Find donors", search)


def _atom_name_element(atom_name):
  atom_name = str(atom_name).strip().upper()
  if not atom_name:
    return ""
  if atom_name[0].isdigit():
    atom_name = atom_name[1:]
  if atom_name.startswith("CL"):
    return "CL"
  if atom_name.startswith("BR"):
    return "BR"
  return atom_name[0]


def _click_spec_to_link_spec(click_spec):
  return [
    _click_spec_chain_id(click_spec),
    _click_spec_res_no(click_spec),
    _click_spec_ins_code(click_spec),
    _click_spec_atom_name(click_spec),
    _click_spec_alt_conf(click_spec),
  ]


def _click_spec_residue_name(click_spec):
  imol = _click_spec_imol(click_spec)
  ch_id = _click_spec_chain_id(click_spec)
  resno = _click_spec_res_no(click_spec)
  ins_code = _click_spec_ins_code(click_spec)
  if imol == -1 or ch_id is False or resno is False:
    return ""
  try:
    return residue_name(imol, ch_id, resno, ins_code)
  except Exception:
    return ""


def _click_spec_summary(click_spec):
  resname = _click_spec_residue_name(click_spec) or "?"
  ch_id = _click_spec_chain_id(click_spec)
  resno = _click_spec_res_no(click_spec)
  atom_name = _click_spec_atom_name(click_spec).strip() or "?"
  return f"{resname} {ch_id}{resno} {atom_name}"


def _click_spec_matches_metal(click_spec, metal_codes):
  return _click_spec_residue_name(click_spec) in metal_codes


def _click_spec_matches_donor(click_spec, donor_elements):
  return _atom_name_element(_click_spec_atom_name(click_spec)) in donor_elements


def _coordination_pair_from_clicks(click_1, click_2, metal_codes, donor_elements):
  first_is_metal = _click_spec_matches_metal(click_1, metal_codes)
  second_is_metal = _click_spec_matches_metal(click_2, metal_codes)
  first_is_donor = _click_spec_matches_donor(click_1, donor_elements)
  second_is_donor = _click_spec_matches_donor(click_2, donor_elements)

  if first_is_metal and second_is_donor:
    return click_1, click_2
  if second_is_metal and first_is_donor:
    return click_2, click_1
  return None, None


def _start_coordination_link_clicks(label, metal_codes, metal_label, donor_elements, donor_label, distance):
  """Start a two-click workflow for ad hoc coordination links."""
  add_status_bar_text(f"{label}: click the metal atom and donor atom")

  def make_link_from_clicks(*args):
    click_1 = args[0]
    click_2 = args[1]
    imol_1 = _click_spec_imol(click_1)
    imol_2 = _click_spec_imol(click_2)
    if imol_1 == -1 or imol_2 == -1 or imol_1 != imol_2:
      info_dialog(f"{label} requires two atoms in the same model molecule.")
      return 0

    metal_click, donor_click = _coordination_pair_from_clicks(click_1, click_2, metal_codes, donor_elements)
    if metal_click is None:
      metal_code_text = "/".join(metal_codes)
      donor_text = "/".join(donor_elements)
      info_dialog(
        f"{label} requires one {metal_label} atom ({metal_code_text}) and one {donor_label} atom ({donor_text}).\n\n"
        f"You clicked:\n{_click_spec_summary(click_1)}\n{_click_spec_summary(click_2)}"
      )
      return 0

    imol = _click_spec_imol(metal_click)
    metal_spec = _click_spec_to_link_spec(metal_click)
    donor_spec = _click_spec_to_link_spec(donor_click)
    coot.make_link_py(imol, metal_spec, donor_spec, "dummy", distance)
    add_status_bar_text(f"{label}: added link ({distance:.3f} A)")
    return 1

  _script_user_defined_click(2, make_link_from_clicks)
  return 1


def make_custom_link():
  """Pick any two atoms in one model, then enter the target link distance."""
  add_status_bar_text("Make custom link: click two atoms")

  def picked(*args):
    if len(args) < 2:
      info_dialog("Make custom link did not receive two atom picks.")
      return 0

    click_1, click_2 = args[0], args[1]
    imol_1 = _click_spec_imol(click_1)
    imol_2 = _click_spec_imol(click_2)
    if imol_1 == -1 or imol_2 == -1 or imol_1 != imol_2:
      info_dialog("Make custom link requires two atoms in the same model molecule.")
      return 0

    atom_1 = _click_spec_to_link_spec(click_1)
    atom_2 = _click_spec_to_link_spec(click_2)
    summary_1 = _click_spec_summary(click_1)
    summary_2 = _click_spec_summary(click_2)

    def submit_distance(value):
      try:
        distance = float(str(value).strip())
      except Exception:
        info_dialog("Custom link distance must be numeric.")
        return 0
      if distance <= 0:
        info_dialog("Custom link distance must be positive.")
        return 0
      try:
        coot.make_link_py(imol_1, atom_1, atom_2, "dummy", distance)
      except Exception as error:
        info_dialog(f"Could not create custom link.\n\n{error}")
        return 0
      add_status_bar_text(
        f"Custom link: {_click_spec_summary(click_1)} - {_click_spec_summary(click_2)} ({distance:.3f} A)"
      )
      return 1

    generic_single_entry(
      f"Custom link target distance (A)\n{summary_1}  -  {summary_2}",
      "2.000",
      "Make link",
      submit_distance,
    )
    return 1

  try:
    _script_user_defined_click(2, picked)
  except Exception as error:
    info_dialog(f"Could not start custom-link atom picking.\n\n{error}")
    return 0
  return 1


def _make_coordination_link(label, metal_codes, metal_label, donor_elements, donor_label, distance, range_text=""):
  """Prompt for a target distance, then start the two-click link helper."""
  def submit_distance(value):
    try:
      chosen_distance = float(str(value).strip())
    except Exception:
      info_dialog(f"{label} requires a numeric distance.")
      return 0
    if chosen_distance <= 0:
      info_dialog(f"{label} requires a positive distance.")
      return 0
    return _start_coordination_link_clicks(
      label,
      metal_codes,
      metal_label,
      donor_elements,
      donor_label,
      chosen_distance,
    )

  generic_single_entry(
    f"{label}: target distance (A){' - ' + range_text if range_text else ''}",
    f"{distance:.3f}",
    "Click two atoms",
    submit_distance,
  )
  return 1


def _parse_phenix_metal_edits(text):
  """Read the explicit bond-add subset of PHIL without a Phenix dependency."""
  import shlex
  scopes, bonds, current = [], [], None
  allowed = {"action", "atom_selection_1", "atom_selection_2", "distance_ideal", "sigma"}
  for line_number, raw in enumerate(text.splitlines(), 1):
    lexer = shlex.shlex(raw, posix=True)
    lexer.whitespace_split = True
    tokens = list(lexer)
    if not tokens:
      continue
    line = raw.strip()
    if len(tokens) == 2 and tokens[1] == "{":
      scope = tokens[0]
      if scope not in ("refinement.geometry_restraints.edits", "refinement", "geometry_restraints", "edits", "bond"):
        raise ValueError(f"Line {line_number}: unsupported scope {scope!r}.")
      if current is not None:
        raise ValueError(f"Line {line_number}: nested bond block.")
      scopes.append(scope)
      if scope == "bond":
        current = {}
    elif tokens == ["}"]:
      if not scopes:
        raise ValueError(f"Line {line_number}: unmatched closing brace.")
      if scopes.pop() == "bond":
        required = allowed - {"action"}
        if not required.issubset(current):
          raise ValueError(f"Bond ending at line {line_number}: missing {sorted(required - current.keys())}.")
        if current.get("action", "add") not in ("add", "*add"):
          raise ValueError(f"Line {line_number}: only action=add is supported.")
        for field in ("distance_ideal", "sigma"):
          value = float(current[field])
          if not math.isfinite(value) or value <= 0:
            raise ValueError(f"Line {line_number}: {field} must be finite and positive.")
          current[field] = value
        bonds.append(current)
        current = None
    else:
      if current is None or "=" not in line:
        raise ValueError(f"Line {line_number}: expected a bond property or scope.")
      key, value = line.split("=", 1)
      key = key.strip()
      if key not in allowed or key in current:
        raise ValueError(f"Line {line_number}: unsupported or repeated property {key!r}.")
      # Preserve quotes in selections (e.g. chain ' ' or name ' CA ').
      current[key] = value.strip()
      if not key.startswith("atom_selection"):
        current[key] = " ".join(shlex.split(value, comments=True))
  if scopes or not bonds:
    raise ValueError("Unclosed scope or no bond edits found.")
  return bonds


def _phenix_metal_atom_spec(imol, selection, residue_index=None):
  """Resolve a restricted Phenix selection to exactly one Coot atom."""
  import shlex
  tokens = shlex.split(selection, comments=True)
  if len(tokens) == 1:
    tokens = shlex.split(tokens[0])
  fields = {}
  while tokens:
    if len(tokens) < 2:
      raise ValueError(f"Incomplete atom selection: {selection}")
    key, value = tokens[:2]
    tokens = tokens[2:]
    if key not in ("name", "chain", "resname", "resseq", "icode", "altloc") or key in fields:
      raise ValueError(f"Unsupported selection term {key!r}: {selection}")
    if any(char in value for char in "*?()"):
      raise ValueError(f"Only explicit atom selections are supported: {selection}")
    fields[key] = value
    if tokens:
      if tokens.pop(0).lower() != "and" or not tokens:
        raise ValueError(f"Only AND selections are supported: {selection}")
  if not {"name", "chain", "resseq"}.issubset(fields):
    raise ValueError(f"Selection requires name, chain and resseq: {selection}")
  chain, resno = fields["chain"], int(fields["resseq"])
  matches = []
  if residue_index is None:
    residue_index = _residue_identity_index(_model_coordinate_snapshot(imol))
  for spec, record in residue_index.get((chain, resno), []):
    ch, rn, ins = spec
    if ch != chain or rn != resno or ("icode" in fields and ins != fields["icode"].strip()):
      continue
    if "resname" in fields and record["name"].strip() != fields["resname"].strip():
      continue
    for atom in record["atoms"]:
      name, alt = atom[0][:2]
      if name.strip() == fields["name"].strip() and ("altloc" not in fields or alt.strip() == fields["altloc"].strip()):
        matches.append([ch, rn, ins, name, alt])
  if len(matches) != 1:
    raise ValueError(f"Selection matched {len(matches)} atoms (expected 1): {selection}")
  return matches[0]


def _residue_identity_index(snapshot):
  index = {}
  for spec, record in snapshot.items():
    index.setdefault(spec[:2], []).append((spec, record))
  return index


def _phenix_link_warning(imol, a, b, snapshot=None):
  """Conservative coordination heuristic, not a full chemical validation."""
  metals = set("LI NA K RB CS FR BE MG CA SR BA RA SC TI V CR MN FE CO NI CU ZN "
               "Y ZR NB MO TC RU RH PD AG CD HF TA W RE OS IR PT AU HG AL GA IN TL "
               "SN PB BI LA CE PR ND PM SM EU GD TB DY HO ER TM YB LU AC TH PA U NP PU AM CM BK CF ES FM MD NO LR".split())
  elements = []
  for spec in (a, b):
    atoms = (snapshot.get(tuple(spec[:3]), {}).get("atoms", []) if snapshot is not None
             else residue_info_py(imol, *spec[:3]) or [])
    matches = [atom for atom in atoms
               if list(atom[0][:2]) == spec[3:]]
    if len(matches) != 1:
      raise ValueError(f"Atom is no longer uniquely available: {spec}")
    elements.append(str(matches[0][1][2]).strip().upper())
  first, second = elements
  if first in metals and second in metals:
    return f"Metal-metal link ({first}-{second})"
  if ((first in metals and second in {"N", "O", "S"}) or
      (second in metals and first in {"N", "O", "S"})):
    return None
  return f"Not a usual metal-N/O/S coordination pair ({first or '?'}-{second or '?'})"


def _apply_phenix_metal_links(imol, prepared):
  if not valid_model_molecule_qm(imol):
    info_dialog("Import cancelled: the target model is no longer available.")
    return 0
  # Revalidate after the asynchronous review, before making any changes.
  try:
    snapshot = _model_coordinate_snapshot(imol)
    for a, b, distance in prepared:
      _phenix_link_warning(imol, a, b, snapshot=snapshot)
  except (ValueError, TypeError, IndexError, RuntimeError) as error:
    info_dialog(f"Import cancelled; no links added.\n\n{error}")
    return 0
  added = 0
  try:
    for a, b, distance in prepared:
      coot.make_link_py(imol, a, b, "dummy", distance)
      added += 1
  except Exception as error:
    info_dialog(f"Import stopped after {added} link calls; earlier changes remain in the model.\n\n{error}")
    return added
  message = f"Submitted {added} Phenix links to Coot for molecule #{imol}."
  print(message)
  add_status_bar_text(message)
  info_dialog(message + "\n\nAtom specifications and ideal distances were passed to make_link. "
              "Sigmas are not imported, and no extra restraints were added. "
              "Coot may use its built-in metal-link distances during refinement.")
  return added


def _review_phenix_metal_links(imol, prepared, warnings_to_review, discarded=None):
  discarded = set() if discarded is None else discarded
  if not warnings_to_review:
    return _apply_phenix_metal_links(
      imol, [link for index, link in enumerate(prepared) if index not in discarded])
  index, warning = warnings_to_review[0]
  a, b, distance = prepared[index]
  def describe(spec):
    chain, number, ins, name, alt = spec
    return f"{chain}:{number}{ins} {name.strip()}" + (f" (alt {alt})" if alt else "")
  def proceed(discard):
    remaining_discarded = discarded | {index} if discard else discarded
    _review_phenix_metal_links(imol, prepared, warnings_to_review[1:], remaining_discarded)
  generic_confirm_dialog(
    "Review unusual Phenix link",
    f"Molecule #{imol}: {warning}\n\n{describe(a)} -- {describe(b)}\n"
    f"Ideal distance: {distance:g} A\n\n"
    "Keep this link or discard it? No links are created until review finishes. "
    "Close this window to cancel the entire import.",
    "Discard link", lambda: proceed(True), "Keep link", lambda: proceed(False))
  return None


def import_phenix_metal_edits(imol, file_name):
  """Create Coot links from Phenix atom selections and ideal distances."""
  try:
    imol = int(imol)
    if not valid_model_molecule_qm(imol):
      raise ValueError("Choose a valid model molecule number.")
    with open(os.path.expanduser(file_name.strip()), encoding="utf-8") as handle:
      bonds = _parse_phenix_metal_edits(handle.read())
    prepared, seen = [], {}
    snapshot = _model_coordinate_snapshot(imol)
    residue_index = _residue_identity_index(snapshot)
    for bond in bonds:
      a = _phenix_metal_atom_spec(imol, bond["atom_selection_1"], residue_index=residue_index)
      b = _phenix_metal_atom_spec(imol, bond["atom_selection_2"], residue_index=residue_index)
      if a == b:
        raise ValueError("A bond cannot join an atom to itself.")
      pair = tuple(sorted((tuple(a), tuple(b))))
      values = bond["distance_ideal"]
      if pair in seen:
        if seen[pair] != values:
          raise ValueError(f"Conflicting edits for atom pair {pair}.")
        continue
      seen[pair] = values
      prepared.append((a, b, values))
    warnings_to_review = []
    for index, (a, b, distance) in enumerate(prepared):
      warning = _phenix_link_warning(imol, a, b, snapshot=snapshot)
      if warning:
        warnings_to_review.append((index, warning))
  except (OSError, ValueError, TypeError, IndexError, RuntimeError) as error:
    info_dialog(f"Phenix metal edits import cancelled; no links added.\n\n{error}")
    return 0
  return _review_phenix_metal_links(imol, prepared, warnings_to_review)


def _links_containing_residue(imol, residue_spec):
  """Read LINK endpoints without counting ordinary dictionary-defined bonds."""
  spec = tuple(residue_spec)
  return [link for link in coot.link_info_py(imol) or []
          if any(tuple(atom[-5:-2]) == spec for atom in link[1:3])]


def remove_links_for_active_residue():
  residue = _active_residue_or_status()
  if not residue:
    return None
  imol, chain, number, ins = residue[:4]
  spec = [chain, number, ins]
  if not all(callable(getattr(coot, name, None)) for name in (
      "link_info_py", "delete_links_containing_residue_py")):
    info_dialog("This Coot build cannot inspect and remove residue links.")
    return None
  try:
    if coot.n_models(imol) != 1:
      raise ValueError("Please split multi-model structures before removing residue links.")
    name = residue_name(imol, *spec)
    links = _links_containing_residue(imol, spec)
  except Exception as error:
    info_dialog(f"Could not inspect residue links.\n\n{error}")
    return None
  label = f"{chain or '(blank)'}:{number}{ins} {name}"
  if not links:
    add_status_bar_text(f"No LINK records involve {label}")
    return None

  def remove_confirmed():
    try:
      if (not coot.is_valid_model_molecule(imol) or coot.n_models(imol) != 1 or
          residue_name(imol, *spec) != name or _links_containing_residue(imol, spec) != links):
        info_dialog("The model or residue links changed. Please reopen Remove links to review them.")
        return None
      with _grouped_model_edit(imol):
        remaining = links
        while remaining:
          # Native deletion can skip links when its loop compacts the LINK list.
          coot.delete_links_containing_residue_py(imol, spec)
          after = _links_containing_residue(imol, spec)
          if len(after) >= len(remaining):
            raise RuntimeError("Coot did not remove all links. Some links may have been removed; use Coot Undo if needed.")
          remaining = after
    except Exception as error:
      info_dialog(f"Could not finish removing residue links.\n\n{error}")
      return None
    add_status_bar_text(f"Removed {len(links)} LINK record(s) involving {label}")
    return len(links)

  return generic_confirm_dialog(
    "Remove links for active residue",
    f"{_molecule_id_name_label(imol)}\n{label}\n\n"
    f"Remove all {len(links)} LINK record(s) involving this residue?\n"
    "This includes coordination and covalent links. Atoms and ordinary residue bonds are unchanged.",
    "Cancel", lambda: None, "Remove links", remove_confirmed,
  )


def import_phenix_metal_edits_dialog():
  imol = _active_molecule_or_status()
  if imol is None:
    return
  initial_directory = None
  model_path = molecule_name(imol)
  if isinstance(model_path, str) and model_path:
    directory = os.path.dirname(os.path.abspath(os.path.expanduser(model_path)))
    if os.path.isdir(directory):
      initial_directory = directory
  generic_single_entry_with_file_browse(
    f"Import Phenix metal edits into active molecule #{imol}", "",
    "Choose Phenix metal edits", "Import",
    lambda file_name: import_phenix_metal_edits(imol, file_name),
    save_file=False,
    initial_directory=initial_directory,
    all_files=True,
  )


def add_coordination_link_menu_entries(menu, entries):
  """Recursively populate the Coordination links submenu."""
  if menu is None:
    return None
  for label, value in entries:
    if isinstance(value, list):
      submenu = Gio.Menu.new()
      menu.append_submenu(label, submenu)
      add_coordination_link_menu_entries(submenu, value)
    else:
      add_simple_coot_menu_menuitem(
        menu,
        label,
        lambda func, params=value, item_label=label: _make_coordination_link(
          item_label,
          params["metal_codes"],
          params["metal_label"],
          params["donor_elements"],
          params["donor_label"],
          params["distance"],
          params.get("range_text", ""),
        ),
      )


def _replace_active_residue_with_monomer(required_resname, required_label, target_resname, modification_label):
  """Apply a replace-residue style covalent modification to the active residue."""
  residue = active_residue()
  if not residue:
    info_dialog(f"{modification_label} requires an active residue.")
    return 0

  mol_id, ch_id, resno, ins_code = residue[:4]
  current_resname = residue_name(mol_id, ch_id, resno, ins_code)
  if current_resname != required_resname:
    info_dialog(f"{modification_label} requires the active residue to be {required_label}.")
    return 0

  # Coot's built-in "Replace residue" GUI path uses mutate_by_overlap(),
  # which supports nonstandard peptide-like residue monomers such as P1L/LLP.
  if ins_code:
    info_dialog(f"{modification_label} is not currently supported for residues with insertion codes.")
    return 0

  status = coot.mutate_by_overlap(mol_id, ch_id, resno, target_resname)
  if status != 1:
    info_dialog(f"Failed to apply {modification_label} to the active {required_label} residue.")
    return 0
  return 1


def _coot_py_residue_spec_to_tuple(res_spec_py):
  return _coot_residue_spec_from_spec(res_spec_py)


def _regularize_linked_residues(imol, base_spec, new_res_spec_py):
  new_spec = _coot_py_residue_spec_to_tuple(new_res_spec_py)
  if not new_spec:
    return 0

  regularize = getattr(coot, "regularize_residues_py", None)
  if not callable(regularize):
    add_status_bar_text("This Coot build cannot regularize an explicit linked residue pair.")
    return 0
  status = regularize(imol, [list(base_spec), list(new_spec)])
  if status:
    accept_regularizement()
  return status


def _atom_xyz_map(imol, chain_id, resno, ins_code):
  atom_info = residue_info_py(imol, chain_id, resno, ins_code)
  if not isinstance(atom_info, list):
    return {}

  coords = {}
  for atom in atom_info:
    try:
      atom_name = atom[0][0]
      xyz = atom[2]
      coords[atom_name] = xyz
    except Exception:
      pass
  return coords


def _first_residue_spec(imol):
  for chain_id in chain_ids(imol):
    resno = seqnum_from_serial_number(imol, chain_id, 0)
    if resno != -10000:
      return (chain_id, resno, "")
  return None


def _find_unused_resno_near(imol, chain_id, start_resno):
  resno = start_resno
  while resno < start_resno + 1000:
    try:
      if not residue_name(imol, chain_id, resno, ""):
        return resno
    except Exception:
      return resno
    resno += 1
  return None


def _merge_ligand_and_make_named_link(required_resname, required_label, ligand_comp_id, link_name,
                                      ligand_atom_name, residue_atom_name, bond_length, modification_label):
  """Import a ligand monomer, merge it into the active model, then create a named link.

  This is used for modifications that are represented by a separate linked
  component rather than a single replacement residue.
  """
  residue = active_residue()
  if not residue:
    info_dialog(f"{modification_label} requires an active residue.")
    return 0

  mol_id, ch_id, resno, ins_code = residue[:4]
  current_resname = residue_name(mol_id, ch_id, resno, ins_code)
  if current_resname != required_resname:
    info_dialog(f"{modification_label} requires the active residue to be {required_label}.")
    return 0

  residue_centre = residue_centre_py(mol_id, ch_id, resno, ins_code)
  if not isinstance(residue_centre, (list, tuple)) or len(residue_centre) != 3:
    info_dialog(f"Failed to locate the active {required_label} residue for {modification_label}.")
    return 0

  ligand_imol = coot.get_monomer(ligand_comp_id)
  if not valid_model_molecule_qm(ligand_imol):
    info_dialog(f"Failed to import {ligand_comp_id} for {modification_label}.")
    return 0

  try:
    ligand_spec = _first_residue_spec(ligand_imol)
    if not ligand_spec:
      info_dialog(f"Failed to prepare {ligand_comp_id} for {modification_label}.")
      return 0

    ligand_chain_id, ligand_resno, ligand_ins_code = ligand_spec
    ligand_coords = _atom_xyz_map(ligand_imol, ligand_chain_id, ligand_resno, ligand_ins_code)
    if not ligand_coords:
      info_dialog(f"Failed to read {ligand_comp_id} coordinates for {modification_label}.")
      return 0

    xs = [xyz[0] for xyz in ligand_coords.values()]
    ys = [xyz[1] for xyz in ligand_coords.values()]
    zs = [xyz[2] for xyz in ligand_coords.values()]
    ligand_centre = [sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs)]
    coot.translate_molecule_by(ligand_imol,
                               residue_centre[0] - ligand_centre[0],
                               residue_centre[1] - ligand_centre[1],
                               residue_centre[2] - ligand_centre[2])

    ligand_resno_new = _find_unused_resno_near(mol_id, ch_id, resno + 1)
    if ligand_resno_new is None:
      info_dialog(f"Failed to find a residue number for {ligand_comp_id} in chain {ch_id}.")
      return 0

    set_merge_molecules_ligand_spec([ch_id, ligand_resno_new, ""])
    merge_result = merge_molecules([ligand_imol], mol_id)
    if not isinstance(merge_result, (list, tuple)) or not merge_result or merge_result[0] != 1:
      info_dialog(f"Failed to merge {ligand_comp_id} for {modification_label}.")
      return 0

    residue_spec = [ch_id, resno, ins_code, residue_atom_name, ""]
    ligand_spec = [ch_id, ligand_resno_new, "", ligand_atom_name, ""]
    coot.make_link_py(mol_id, residue_spec, ligand_spec, link_name, bond_length)
    _regularize_linked_residues(mol_id, (ch_id, resno, ins_code), [ch_id, ligand_resno_new, ""])
    return 1
  finally:
    if valid_model_molecule_qm(ligand_imol):
      coot.close_molecule(ligand_imol)


def _merge_ligand_and_make_named_link_for_active_residue(link_specs_by_resname, modification_label):
  """Dispatch a merge-link modification based on the active residue type."""
  residue = active_residue()
  if not residue:
    info_dialog(f"{modification_label} requires an active residue.")
    return 0

  mol_id, ch_id, resno, ins_code = residue[:4]
  current_resname = residue_name(mol_id, ch_id, resno, ins_code)
  link_spec = link_specs_by_resname.get(current_resname)
  if not link_spec:
    valid_labels = [spec["required_label"] for _resname, spec in sorted(link_specs_by_resname.items())]
    valid_text = " or ".join(valid_labels)
    info_dialog(f"{modification_label} requires the active residue to be {valid_text}.")
    return 0

  return _merge_ligand_and_make_named_link(
    link_spec["required_resname"],
    link_spec["required_label"],
    link_spec["ligand_comp_id"],
    link_spec["link_name"],
    link_spec["ligand_atom_name"],
    link_spec["residue_atom_name"],
    link_spec["bond_length"],
    modification_label,
  )


COVALENT_MODIFICATION_MENU = [
  ("Palmitoylation (Cys)", {
    "mode": "replace",
    "required_resname": "CYS",
    "required_label": "Cys",
    "target_resname": "P1L",
    "modification_label": "Palmitoylation",
  }),
  ("BME adduct (Cys)", {
    "mode": "merge_link",
    "required_resname": "CYS",
    "required_label": "Cys",
    "ligand_comp_id": "BME",
    "link_name": "CYS-BME",
    "ligand_atom_name": " S2 ",
    "residue_atom_name": " SG ",
    "bond_length": 2.023,
    "modification_label": "BME adduct",
  }),
  ("PLP linkage (Lys)", {
    "mode": "merge_link",
    "required_resname": "LYS",
    "required_label": "Lys",
    "ligand_comp_id": "PLP",
    "link_name": "LYS-PLP",
    "ligand_atom_name": " C4A",
    "residue_atom_name": " NZ ",
    "bond_length": 1.270,
    "modification_label": "PLP linkage",
  }),
  ("Monomethyl-Lys (Lys)", {
    "mode": "replace",
    "required_resname": "LYS",
    "required_label": "Lys",
    "target_resname": "MLZ",
    "modification_label": "Monomethyl-Lys",
  }),
  ("Dimethyl-Lys (Lys)", {
    "mode": "replace",
    "required_resname": "LYS",
    "required_label": "Lys",
    "target_resname": "MLY",
    "modification_label": "Dimethyl-Lys",
  }),
  ("Trimethyl-Lys (Lys)", {
    "mode": "replace",
    "required_resname": "LYS",
    "required_label": "Lys",
    "target_resname": "M3L",
    "modification_label": "Trimethyl-Lys",
  }),
  ("Acetyl-Lys (Lys)", {
    "mode": "replace",
    "required_resname": "LYS",
    "required_label": "Lys",
    "target_resname": "ALY",
    "modification_label": "Acetyl-Lys",
  }),
  ("Carboxy-Lys (Lys)", {
    "mode": "replace",
    "required_resname": "LYS",
    "required_label": "Lys",
    "target_resname": "KCX",
    "modification_label": "Carboxy-Lys",
  }),
  ("Retinal linkage (Lys)", {
    "mode": "merge_link",
    "required_resname": "LYS",
    "required_label": "Lys",
    "ligand_comp_id": "RET",
    "link_name": "LYS-RET",
    "ligand_atom_name": " C15",
    "residue_atom_name": " NZ ",
    "bond_length": 1.267,
    "modification_label": "Retinal linkage",
  }),
  ("Phosphoserine (Ser)", {
    "mode": "replace",
    "required_resname": "SER",
    "required_label": "Ser",
    "target_resname": "SEP",
    "modification_label": "Phosphoserine",
  }),
  ("Phosphothreonine (Thr)", {
    "mode": "replace",
    "required_resname": "THR",
    "required_label": "Thr",
    "target_resname": "TPO",
    "modification_label": "Phosphothreonine",
  }),
  ("O-linked NAG (Ser/Thr)", {
    "mode": "merge_link_by_resname",
    "modification_label": "O-linked NAG",
    "link_specs_by_resname": {
      "SER": {
        "required_resname": "SER",
        "required_label": "Ser",
        "ligand_comp_id": "NAG",
        "link_name": "pyr-SER",
        "ligand_atom_name": " C1 ",
        "residue_atom_name": " OG ",
        "bond_length": 1.404,
      },
      "THR": {
        "required_resname": "THR",
        "required_label": "Thr",
        "ligand_comp_id": "NAG",
        "link_name": "pyr-THR",
        "ligand_atom_name": " C1 ",
        "residue_atom_name": " OG1",
        "bond_length": 1.391,
      },
    },
  }),
  ("Phosphotyrosine (Tyr)", {
    "mode": "replace",
    "required_resname": "TYR",
    "required_label": "Tyr",
    "target_resname": "PTR",
    "modification_label": "Phosphotyrosine",
  }),
  ("Selenomethionine (Met)", {
    "mode": "replace",
    "required_resname": "MET",
    "required_label": "Met",
    "target_resname": "MSE",
    "modification_label": "Selenomethionine",
  }),
  ("Methionine sulfoxide (Met)", {
    "mode": "replace",
    "required_resname": "MET",
    "required_label": "Met",
    "target_resname": "SME",
    "modification_label": "Methionine sulfoxide",
  }),
  ("Hydroxyproline (Pro)", {
    "mode": "replace",
    "required_resname": "PRO",
    "required_label": "Pro",
    "target_resname": "HYP",
    "modification_label": "Hydroxyproline",
  }),
  ("Citrullination (Arg)", {
    "mode": "replace",
    "required_resname": "ARG",
    "required_label": "Arg",
    "target_resname": "CIR",
    "modification_label": "Citrullination",
  }),
  ("Methyl-Arg (Arg)", {
    "mode": "replace",
    "required_resname": "ARG",
    "required_label": "Arg",
    "target_resname": "AGM",
    "modification_label": "Methyl-Arg",
  }),
  ("S-hydroxycysteine (Cys)", {
    "mode": "replace",
    "required_resname": "CYS",
    "required_label": "Cys",
    "target_resname": "CSO",
    "modification_label": "S-hydroxycysteine",
  }),
]


def _apply_covalent_modification(spec):
  """Dispatch a covalent-modification menu entry from a single data record."""
  if spec["mode"] == "replace":
    return _replace_active_residue_with_monomer(
      spec["required_resname"],
      spec["required_label"],
      spec["target_resname"],
      spec["modification_label"],
    )
  if spec["mode"] == "merge_link":
    return _merge_ligand_and_make_named_link(
      spec["required_resname"],
      spec["required_label"],
      spec["ligand_comp_id"],
      spec["link_name"],
      spec["ligand_atom_name"],
      spec["residue_atom_name"],
      spec["bond_length"],
      spec["modification_label"],
    )
  if spec["mode"] == "merge_link_by_resname":
    return _merge_ligand_and_make_named_link_for_active_residue(
      spec["link_specs_by_resname"],
      spec["modification_label"],
    )
  raise ValueError("Unknown covalent modification mode: {0}".format(spec["mode"]))


def add_covalent_modification_menu_entries(menu):
  """Populate the Build -> Covalent modifications submenu."""
  if menu is None:
    return None

  for label, spec in COVALENT_MODIFICATION_MENU:
    add_simple_coot_menu_menuitem(
      menu,
      label,
      lambda func, item_spec=spec: _apply_covalent_modification(item_spec),
    )


#Colors subset of protein residues red, provided by user as string of single-letter ids.
def color_protein_residue_subset():
  residue = _active_residue_or_status()
  if not residue:
    return None
  mol_id = residue[0]
  def color_from_string(X):
    entry=str(X).upper() #capitalize
    aa_dic={'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS','E':'GLU','Q':'GLN','G':'GLY','H':'HIS','I':'ILE','L':'LEU','K':'LYS','M':'MET','F':'PHE','P':'PRO','S':'SER','T':'THR','W':'TRP','Y':'TYR','V':'VAL','X':'UNK'}
    highlight_colour=34 #red
    blank_colour=0 #light grey
    residue_list=[]
    resname_list=[]
    for char in entry:
      if (aa_dic.get(char,0)!=0): #make a list of 3-letter resnames from user supplied string by checking aa_dic
        resname_list.append(aa_dic.get(char,0))
    for ch_id in chain_ids(mol_id):
      for sn in _iter_chain_serial_numbers(mol_id, ch_id):
        resname_here=resname_from_serial_number(mol_id,ch_id,sn)
        for resname in resname_list:
          if resname==resname_here:
            resn=seqnum_from_serial_number(mol_id,ch_id,sn)
            ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
            residue_list.append(([ch_id,resn,ins_id],highlight_colour))
        if resname_here not in resname_list:
          resn=seqnum_from_serial_number(mol_id,ch_id,sn)
          ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
          residue_list.append(([ch_id,resn,ins_id],blank_colour))
    _apply_user_defined_residue_colours(mol_id, [], residue_list)
  generic_single_entry("Residues to color? Single letter code (e.g. DE or de will color Asp/Glu)","A","Color entered residue types!",color_from_string)

def _place_new_secondary_structure(prompt_label, button_label, phi, psi):
  """Create a new poly-Ala secondary-structure fragment with idealized torsions."""
  def place_from_entry(n):
    count = _parse_int_entry(n, f"{button_label} requires a positive integer.", minimum_value=1)
    if count is None:
      return 0
    mol_id = get_monomer_no_H("ALA")
    if not valid_model_molecule_qm(mol_id):
      info_dialog(f"Failed to create the starting monomer for {button_label}.")
      return 0
    ch_id=chain_ids(mol_id)[0]
    res_no=1
    ins_code=""
    altloc=""
    delete_atom(mol_id,ch_id,1,ins_code," OXT",altloc)
    for _i in range(1, count):
      add_terminal_residue_using_phi_psi(mol_id,ch_id,res_no,"auto",phi,psi)
      res_no=res_no+1
    return 1
  generic_single_entry(prompt_label, "10", button_label, place_from_entry)

#Mutate active chain to entered sequence
default_seq="MAAAA"
def mutate_by_resnum():
  chain_context = _active_chain_context_or_status()
  if not chain_context:
    return None
  mol_id = chain_context["mol_id"]
  ch_id = chain_context["chain_id"]
  def enter_seq(seq):
    global default_seq
    seq=str(seq).upper()
    seq=seq.replace(" ", "")
    seq_dic={}
    len_seq=len(seq)
    n=0
    nmax=len_seq+1
    aa_dic={'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS','E':'GLU',
    'Q':'GLN','G':'GLY','H':'HIS','I':'ILE','L':'LEU','K':'LYS','M':'MET',
    'F':'PHE','P':'PRO','S':'SER','T':'THR','W':'TRP','Y':'TYR','V':'VAL'}
    valid_polymer_resnames = set(aa_dic.values()) | {"MSE"}
    clean_seq_chars=[]
    while (n<len_seq): 
      if ((seq[n].isalpha() and (seq[n] in aa_dic))):
        clean_seq_chars.append(seq[n])
      n=n+1
    seq="".join(clean_seq_chars)
    default_seq=seq
    len_seq=len(seq)
    n=0
    while (n<len_seq):
      value=aa_dic[seq[n]]
      seq_dic[n+1]=value
      n=n+1
    sn=0
    last_sn=chain_n_residues(ch_id,mol_id)-1
    with _grouped_model_edit(mol_id):
      while (sn<=last_sn):
        res=resname_from_serial_number(mol_id,ch_id,sn)
        seqnum=seqnum_from_serial_number(mol_id,ch_id,sn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
        if ((res!=seq_dic.get(seqnum)) and (res in valid_polymer_resnames)):
          if seq_dic.get(seqnum):
            mutate(mol_id,ch_id,seqnum,ins_id,seq_dic.get(seqnum))
            if seq_dic.get(seqnum)!="PRO":
              delete_residue_sidechain(mol_id,ch_id,seqnum,ins_id,0)
        sn=sn+1
  generic_single_entry("Enter raw amino acid sequence (must be complete!)",
  default_seq,"Mutate active chain to match sequence using PDB numbering", enter_seq)

#Make new helix (don't fit)
def place_new_helix():
  return _place_new_secondary_structure(
    "How many residues for helix?",
    "Place helix",
    -57.82,
    -47,
  )
  
#Make new strand (don't fit)
def place_new_strand():
  return _place_new_secondary_structure(
    "How many residues for strand?",
    "Place strand",
    -139,
    135,
  )

#Make new 3-10 helix (don't fit)
def place_new_3_10_helix():
  return _place_new_secondary_structure(
    "How many residues for 3-10 helix?",
    "Place 3-10 helix",
    -49,
    -26,
  )

#Renumber active segment by active residue
def renumber_seg_by_active_res():
  context = _active_segment_context_or_status()
  if not context:
    return None
  current_num = context["resno"]
  mol_id = context["mol_id"]
  ch_id = context["chain_id"]
  selected = set(context["segment_specs"])
  def renum_seg(new_num):
    new_num = _parse_int_entry(new_num, "Renumbering requires an integer residue number.")
    if new_num is None:
      return 0
    specs = [tuple(spec) for spec in _tq_chain_residue_specs(mol_id, ch_id)]
    start, end = min(spec[1] for spec in selected), max(spec[1] for spec in selected)
    if {spec for spec in specs if start <= spec[1] <= end} != selected:
      info_dialog("The numbered range includes residues outside this segment; nothing was renumbered.")
      return 0
    offset = new_num - current_num
    occupied = set(specs) - selected
    if any((chain, number + offset, ins) in occupied for chain, number, ins in selected):
      info_dialog("Renumbering would overlap existing residues; nothing was renumbered.")
      return 0
    if not _renumber_residue_range_or_dialog(mol_id, ch_id, start, end, offset):
      return 0
    _refresh_extra_restraints_display(mol_id)
    return 1
  generic_single_entry("New number for this residue?",
  str(current_num),"Renumber",renum_seg)

def find_sequence_with_entry():
  generic_single_entry("Enter sequence fragment to find\n(X=wildcard, (A/S)=either residue)",
  "MAAAA","Find sequence in active chain",find_sequence_in_current_chain)

def user_defined_add_arbitrary_length_bond_restraint(bond_length=2.0):
  def make_restr(text_list, continue_qm):
    s = "Now click on 2 atoms to define the additional bond restraint"
    add_status_bar_text(s)
    dist = text_list[0]
    try:
      bl = float(dist)
    except (TypeError, ValueError):
      bl = False
      add_status_bar_text("Must define a number for the bond length")
    if bl:
      def make_restr_dist(*args):
        atom_spec_1 = args[0]
        atom_spec_2 = args[1]
        imol = _click_spec_imol(atom_spec_1)
        if imol < 0 or imol != _click_spec_imol(atom_spec_2):
          info_dialog("Extra distance restraints require two atoms in the same model.")
          return
        fields = (_click_spec_chain_id, _click_spec_res_no, _click_spec_ins_code,
                  _click_spec_atom_name, _click_spec_alt_conf)
        add_extra_bond_restraint(imol, *(f(spec) for spec in (atom_spec_1, atom_spec_2)
                                       for f in fields), bl, 0.035)
        if continue_qm:
          user_defined_add_arbitrary_length_bond_restraint(bl)
      _script_user_defined_click(2, make_restr_dist)
  def stay_open(*args):
    pass
  #generic_single_entry("Add a User-defined extra distance restraint",
  #                     "2.0",
  #                     "OK...",
  #                     lambda text: make_restr(text))
  generic_multiple_entries_with_check_button(
    [["Add a User-defined extra distance restraint",
      str(bond_length)]],
    ["Stay open?", lambda active_state: stay_open(active_state)],
    "OK...",
    lambda text, stay_open_qm: make_restr(text, stay_open_qm))


def toggle_refinement_ramachandran_overlay():
  """Persist plot visibility without changing Ramachandran restraints."""
  # No native getter exists; track our changes from Coot's enabled default.
  enabled = not REFINEMENT_RAMACHANDRAN_OVERLAY_ENABLED
  if not _set_refinement_ramachandran_overlay(enabled):
    return None
  saved = _save_display_default("refinement_ramachandran_overlay", enabled)
  add_status_bar_text("Refinement Ramachandran overlay: " + ("on" if enabled else "off") +
                      (" (saved for future sessions)" if saved else " (this session only; could not save)"))
  return enabled


def _build_custom_display_menu(submenu_display):
  """Populate the Display submenu in normal Python rather than exec'd text."""
  add_simple_coot_menu_menuitem(
    submenu_display, "Trackpad controls...",
    lambda func: custom_trackpad_controls_dialog(),
  )
  add_simple_coot_menu_menuitem(
    submenu_display,
    "Default model to high contrast display",
    lambda func: default_models_to_high_contrast_display(),
  )
  add_simple_coot_menu_menuitem(
    submenu_display,
    "Default to ambient maps",
    lambda func: default_maps_to_ambient_display(),
  )
  add_simple_coot_menu_menuitem(
    submenu_display,
    "Default to black background",
    lambda func: default_background_to_black(),
  )
  add_simple_coot_menu_menuitem(
    submenu_display,
    "Reset display defaults (restart required)",
    lambda func: reset_map_model_defaults(),
  )
  add_simple_coot_menu_menuitem(
    submenu_display,
    "All Molecules use \"C-alpha\" Symmetry",
    lambda func: [valid_model_molecule_qm(imol) and symmetry_as_calphas(imol, 1) for imol in molecule_number_list()],
  )
  add_simple_coot_menu_menuitem(
    submenu_display,
    "Toggle Symmetry",
    lambda func: set_show_symmetry_master(not get_show_symmetry()),
  )
  add_simple_coot_menu_menuitem(
    submenu_display,
    "Clear labels and distances",
    lambda func: clear_distances_and_labels(),
  )
  add_simple_coot_menu_menuitem(
    submenu_display,
    "Toggle high-contrast model lighting",
    lambda func: toggle_high_contrast_mode(),
  )
  add_simple_coot_menu_menuitem(
    submenu_display,
    "Toggle Ramachandran overlay during refinement",
    lambda func: toggle_refinement_ramachandran_overlay(),
  )
  add_simple_coot_menu_menuitem(
    submenu_display,
    "Toggle Rama/rotamer spheres during refinement",
    lambda func: toggle_refinement_validation_spheres(),
  )
  add_simple_coot_menu_menuitem(
    submenu_display,
    "Yellow-ify carbons",
    lambda func: yellowify_carbons_in_active_molecule(),
  )
  add_simple_coot_menu_menuitem(
    submenu_display,
    "Switch all mols to CA representation",
    lambda func: all_mols_to_ca(),
  )
  add_simple_coot_menu_menuitem(
    submenu_display,
    "Resample/restyle current map",
    lambda func: resample_active_map_for_em_half_angstrom(),
  )
  add_simple_coot_menu_menuitem(
    submenu_display,
    "Map brightness...",
    lambda func: show_active_map_brightness_slider(),
  )
  add_simple_coot_menu_menuitem(
    submenu_display,
    "Map colour...",
    lambda func: show_active_map_colour_slider(),
  )
  add_simple_coot_menu_menuitem(
    submenu_display,
    "Find sequence in active chain",
    lambda func: find_sequence_with_entry(),
  )
  add_simple_coot_menu_menuitem(
    submenu_display,
    "Sequence context",
    lambda func: sequence_context(),
  )
  add_simple_coot_menu_menuitem(
    submenu_display,
    "Odd residues",
    lambda func: find_odd_residues(),
  )
  add_simple_coot_menu_menuitem(
    submenu_display,
    "Residue annotations...",
    lambda func: residue_annotations_gui(),
  )
  add_simple_coot_menu_menuitem(
    submenu_display,
    "Nearby residue annotations...",
    lambda func: nearby_residue_annotations_gui(),
  )


def _build_custom_user_colour_menu(submenu_colour):
  add_simple_coot_menu_menuitem(
    submenu_colour,
    "Color active mol by rotamer prob (outliers magenta) and missing atoms (blue)",
    lambda func: color_rotamer_outliers_and_missing_atoms_for_active_molecule(),
  )
  add_simple_coot_menu_menuitem(
    submenu_colour,
    "Color active mol by hydrophobics (orange), polars (blue), glys (magenta) and pros (green)",
    lambda func: color_polars_and_hphobs_for_active_molecule(),
  )
  add_simple_coot_menu_menuitem(
    submenu_colour,
    "Color active mol by charge (+ve blue, -ve red)",
    lambda func: color_by_charge_for_active_molecule(),
  )
  add_simple_coot_menu_menuitem(
    submenu_colour,
    "Uncolor other chains in active mol",
    lambda func: uncolor_other_chains(),
  )
  add_simple_coot_menu_menuitem(
    submenu_colour,
    "Color active chain",
    lambda func: color_active_chain(),
  )
  add_simple_coot_menu_menuitem(
    submenu_colour,
    "Color active segment",
    lambda func: colour_active_segment(),
  )
  add_simple_coot_menu_menuitem(
    submenu_colour,
    "Color by protein/nucleic acid",
    lambda func: color_protein_na_for_active_molecule(),
  )
  add_simple_coot_menu_menuitem(
    submenu_colour,
    "Color waters",
    lambda func: color_waters_for_active_molecule(),
  )
  add_simple_coot_menu_menuitem(
    submenu_colour,
    "Colour entered subset of protein residues for active mol",
    lambda func: color_protein_residue_subset(),
  )
  add_simple_coot_menu_menuitem(
    submenu_colour,
    "Color active mol by Ramachandran outliers",
    lambda func: color_by_rama_native_for_active_residue(),
  )
  add_simple_coot_menu_menuitem(
    submenu_colour,
    "Color active mol by density fit",
    lambda func: color_by_density_fit_native_for_active_residue(),
  )
  add_simple_coot_menu_menuitem(
    submenu_colour,
    "Color active mol by NCS difference",
    lambda func: color_by_ncs_difference_for_active_residue(),
  )
  add_simple_coot_menu_menuitem(
    submenu_colour,
    "Color active mol by clash score",
    lambda func: color_by_clash_score_for_active_molecule(),
  )
  add_simple_coot_menu_menuitem(
    submenu_colour,
    "Clear custom colours",
    lambda func: clear_custom_colour_representations_for_active_molecule(),
  )


def _build_custom_fit_menu(submenu_fit):
  trackpad_menu = Gio.Menu.new()
  for scope, label in (("residue", "Active residue..."), ("range", "Clicked residue range..."),
                       ("chain", "Active chain..."), ("molecule", "Active molecule...")):
    add_simple_coot_menu_menuitem(trackpad_menu, label,
                                 lambda _item, choice=scope: trackpad_rotate_translate(choice))
  submenu_fit.append_submenu("Trackpad Rotate/Translate", trackpad_menu)
  add_simple_coot_menu_menuitem(submenu_fit, "Cylinder refine (click range)...", lambda func: cylinder_refine())
  add_simple_coot_menu_menuitem(submenu_fit, "Cylinder refine (alpha helix)...", lambda func: cylinder_refine_alpha_helix())
  add_simple_coot_menu_menuitem(submenu_fit, "Cylinder refine (beta strand)...", lambda func: cylinder_refine_beta_strand())
  add_simple_coot_menu_menuitem(
    submenu_fit,
    "Multimap stepped refinement...",
    lambda func: multimap_refinement_gui(),
  )
  add_simple_coot_menu_menuitem(submenu_fit, "Fit all chains to map", lambda func: rigid_fit_all_chains())
  add_simple_coot_menu_menuitem(submenu_fit, "Fit current chain to map", lambda func: rigid_fit_active_chain())
  add_simple_coot_menu_menuitem(submenu_fit, "Global fit to map...", lambda func: global_fit_map())
  add_simple_coot_menu_menuitem(
    submenu_fit,
    "Jiggle-fit current chain to map (Slow!)",
    lambda func: jiggle_fit_active_chain(),
  )
  add_simple_coot_menu_menuitem(
    submenu_fit,
    "Jiggle-fit current chain to B-smoothed map (Slow!)",
    lambda func: jiggle_fit_active_chain_smooth(),
  )
  add_simple_coot_menu_menuitem(
    submenu_fit,
    "Jiggle-fit all chains to map (very slow!)",
    lambda func: jiggle_fit_all_chains(),
  )
  add_simple_coot_menu_menuitem(
    submenu_fit,
    "Jiggle-fit current mol to map (Slow!)",
    lambda func: jiggle_fit_active_mol(),
  )
  add_simple_coot_menu_menuitem(submenu_fit, "Fit all segments", lambda func: rigid_body_fit_segments())
  add_simple_coot_menu_menuitem(submenu_fit, "Fit this segment", lambda func: fit_this_segment())
  add_simple_coot_menu_menuitem(
    submenu_fit,
    "Smart self restrain active mol...",
    lambda func: prompt_generate_smart_local_extra_restraints(),
  )


def _build_custom_renumber_menu(submenu_renumber):
  add_simple_coot_menu_menuitem(
    submenu_renumber,
    "Renumber active chain by first res",
    lambda func: renumber_by_first_res(),
  )
  add_simple_coot_menu_menuitem(
    submenu_renumber,
    "Renumber active chain by last res",
    lambda func: renumber_by_last_res(),
  )
  add_simple_coot_menu_menuitem(
    submenu_renumber,
    "Renumber active chain by current res",
    lambda func: renumber_by_active_res(),
  )
  add_simple_coot_menu_menuitem(
    submenu_renumber,
    "Renumber from N-term to active residue",
    lambda func: renumber_n_term_segment(),
  )
  add_simple_coot_menu_menuitem(
    submenu_renumber,
    "Renumber from active residue to C-term",
    lambda func: renumber_c_term_segment(),
  )
  add_simple_coot_menu_menuitem(
    submenu_renumber,
    "Renumber segment by active res",
    lambda func: renumber_seg_by_active_res(),
  )


def _build_custom_settings_menu(submenu_settings):
  add_simple_coot_menu_menuitem(
    submenu_settings,
    "Auto-scale B-factor coloring for active mol",
    lambda func: autoscale_b_factor(),
  )
  add_simple_coot_menu_menuitem(
    submenu_settings,
    "Set Bfac for new atoms to mean B for active mol",
    lambda func: set_new_atom_b_fac_to_mean(),
  )
  add_simple_coot_menu_menuitem(
    submenu_settings,
    "Set proportional editing radius",
    lambda func: set_proportional_editing_radius(),
  )


def _build_custom_build_menu(
  submenu_build,
  submenu_common_monomers,
  submenu_common_monomer_favorites,
  submenu_coordination_links,
  submenu_covalent_modifications,
):
  global COMMON_MONOMER_FAVORITES_MENU

  add_simple_coot_menu_menuitem(
    submenu_build,
    "Forced addition of terminal residue (click terminus)",
    lambda func: force_add_terminal_residue(),
  )
  add_simple_coot_menu_menuitem(
    submenu_build,
    "Grow helix (click terminus)",
    lambda func: grow_helix(),
  )
  add_simple_coot_menu_menuitem(
    submenu_build,
    "Grow strand (click terminus)",
    lambda func: grow_strand(),
  )
  add_simple_coot_menu_menuitem(
    submenu_build,
    "Grow parallel strand (click terminus)",
    lambda func: grow_parallel_strand(),
  )
  add_simple_coot_menu_menuitem(
    submenu_build,
    "Grow 3-10 helix (click terminus)",
    lambda func: grow_helix_3_10(),
  )
  add_simple_coot_menu_menuitem(
    submenu_build,
    "Shorten loop by one residue",
    lambda func: shorten_loop(),
  )
  add_simple_coot_menu_menuitem(
    submenu_build,
    "Lengthen loop by one residue",
    lambda func: lengthen_loop(),
  )
  add_simple_coot_menu_menuitem(
    submenu_build,
    "Get fractional coordinates of active atom",
    lambda func: get_fract_coords(),
  )

  submenu_build.append_submenu("Common monomers", submenu_common_monomers)
  submenu_common_monomers.append_submenu("Favorites", submenu_common_monomer_favorites)
  COMMON_MONOMER_FAVORITES_MENU = submenu_common_monomer_favorites
  refresh_common_monomer_favorites_menu()
  add_common_monomer_menu_entries(submenu_common_monomers, COMMON_MONOMER_MENU)

  submenu_build.append_submenu("Coordination links", submenu_coordination_links)
  add_coordination_link_menu_entries(submenu_coordination_links, COORDINATION_LINK_MENU)
  add_simple_coot_menu_menuitem(
    submenu_coordination_links, "Make custom link...",
    lambda *_args: make_custom_link(),
  )
  add_simple_coot_menu_menuitem(
    submenu_coordination_links, "Propose metal links for current residue...",
    lambda func: auto_metal_links_current_residue(),
  )
  add_simple_coot_menu_menuitem(
    submenu_coordination_links, "Import Phenix metal edits...",
    lambda *_args: import_phenix_metal_edits_dialog(),
  )
  add_simple_coot_menu_menuitem(
    submenu_coordination_links, "Remove links for active residue...",
    lambda *_args: remove_links_for_active_residue(),
  )

  submenu_build.append_submenu("Covalent modifications", submenu_covalent_modifications)
  add_covalent_modification_menu_entries(submenu_covalent_modifications)

  add_simple_coot_menu_menuitem(
    submenu_build,
    "Make alpha helix of length n",
    lambda func: place_new_helix(),
  )
  add_simple_coot_menu_menuitem(
    submenu_build,
    "Make 3-10 helix of length n",
    lambda func: place_new_3_10_helix(),
  )
  add_simple_coot_menu_menuitem(
    submenu_build,
    "Build polyala loop (click start,end)",
    lambda func: fit_polyala_gui(),
  )


def _build_custom_mutate_menu(submenu_mutate):
  add_simple_coot_menu_menuitem(
    submenu_mutate, "Compare local sequence registers...",
    lambda func: compare_local_sequence_registers(),
  )
  add_simple_coot_menu_menuitem(
    submenu_mutate, "Shift sequence register of active segment...",
    lambda func: shift_register_active_segment(),
  )
  add_simple_coot_menu_menuitem(
    submenu_mutate, "Shift sequence register (click start and end)...",
    lambda func: shift_register_clicked_range(),
  )
  add_simple_coot_menu_menuitem(
    submenu_mutate,
    "Mutate range to UNK (click start and end)",
    lambda func: mutate_residue_range_by_click_a(),
  )
  add_simple_coot_menu_menuitem(
    submenu_mutate,
    "Mutate range to ALA (click start and end)",
    lambda func: mutate_residue_range_by_click_ala_a(),
  )
  add_simple_coot_menu_menuitem(
    submenu_mutate,
    "Mutate all Mets to MSE",
    lambda func: mutate_all_mets_to_mse(),
  )
  add_simple_coot_menu_menuitem(
    submenu_mutate,
    "Mutate all MSEs to Met",
    lambda func: mutate_all_mse_to_met(),
  )
  add_simple_coot_menu_menuitem(
    submenu_mutate,
    "Mutate active chain to template sequence (numbering must match sequence!)",
    lambda func: mutate_by_resnum(),
  )


def _build_custom_modify_menu(submenu_modify):
  add_simple_coot_menu_menuitem(submenu_modify, "Copy current chain", lambda func: copy_active_chain())
  add_simple_coot_menu_menuitem(submenu_modify, "Cut current chain", lambda func: cut_active_chain())
  add_simple_coot_menu_menuitem(submenu_modify, "Copy active segment", lambda func: copy_active_segment())
  add_simple_coot_menu_menuitem(submenu_modify, "Cut active segment", lambda func: cut_active_segment())
  add_simple_coot_menu_menuitem(
    submenu_modify,
    "Smart copy active non-polymer residue",
    lambda func: smart_copy_active_non_polymer_residue(),
  )
  add_simple_coot_menu_menuitem(
    submenu_modify,
    "Smart paste copied non-polymer residue",
    lambda func: smart_paste_copied_non_polymer_residue(),
  )
  add_simple_coot_menu_menuitem(
    submenu_modify,
    "Copy fragment (click start and end)",
    lambda func: copy_frag_by_click(),
  )
  add_simple_coot_menu_menuitem(
    submenu_modify,
    "Cut fragment (click start and end)",
    lambda func: cut_frag_by_click(),
  )
  add_simple_coot_menu_menuitem(
    submenu_modify,
    "Copy active chain to NCS equivs",
    lambda func: copy_ncs_chain_from_active(),
  )
  add_simple_coot_menu_menuitem(
    submenu_modify,
    "Delete active segment",
    lambda func: delete_active_segment(),
  )
  add_simple_coot_menu_menuitem(
    submenu_modify,
    "Merge chains (click two; 2nd into 1st)",
    lambda func: merge_chains(),
  )


def _build_custom_maps_menu(submenu_maps):
  add_simple_coot_menu_menuitem(
    submenu_maps,
    "Go to center of scrollable map",
    lambda func: goto_center_of_map(),
  )
  add_simple_coot_menu_menuitem(
    submenu_maps,
    "Set refinement map to scrollable map",
    lambda func: set_map_to_scrollable_map(),
  )
  add_simple_coot_menu_menuitem(
    submenu_maps,
    "Make masked map (EM)",
    lambda func: make_masked_map_em(),
  )


def _custom_menu_search_items(menu, path=(), ancestors=()):
  """Read live GMenu actions, including dynamic submenus, without copying callbacks."""
  if any(menu == parent for parent in ancestors) or len(ancestors) > 20:
    return []
  items = []
  def attribute(index, name):
    value = menu.get_item_attribute_value(index, name, None)
    return value.unpack() if value is not None else None
  for index in range(menu.get_n_items()):
    label = attribute(index, "label") or ""
    # GMenu uses single underscores for mnemonics and doubles literal ones.
    label = label.replace("__", "\0").replace("_", "").replace("\0", "_")
    submenu = menu.get_item_link(index, "submenu")
    section = menu.get_item_link(index, "section")
    if submenu is not None or section is not None:
      child = submenu if submenu is not None else section
      items.extend(_custom_menu_search_items(child, path + ((label,) if label else ()), ancestors + (menu,)))
    else:
      action = attribute(index, "action")
      if label and action:
        items.append({"label": label, "path": " > ".join(path), "action": action,
                      "target": menu.get_item_attribute_value(index, "target", None)})
  return items


def _custom_menu_search_matches(menu, query):
  terms = query.casefold().split()
  if not terms:
    return []
  items = [item for item in _custom_menu_search_items(menu)
           if all(term in (item["path"] + " " + item["label"]).casefold() for term in terms)]
  phrase = " ".join(terms)
  return sorted(items, key=lambda item: (not item["label"].casefold().startswith(phrase),
                                         item["path"].casefold(), item["label"].casefold()))


def _menu_search_items(custom_menu, native_menu=None):
  """Combine live menus for search without changing the toolbar's Custom catalog."""
  custom = _custom_menu_search_items(custom_menu)
  def identity(item):
    target = item["target"]
    return item["action"], None if target is None else (target.get_type_string(), target.print_(True))
  custom_actions = {identity(item) for item in custom}
  items = [dict(item, source="Custom") for item in custom]
  if native_menu is not None:
    items.extend(dict(item, source="Coot") for item in _custom_menu_search_items(native_menu)
                 if identity(item) not in custom_actions)
  return items


def _menu_search_matches(items, query):
  terms = query.casefold().split()
  if not terms:
    return []
  phrase = " ".join(terms)
  return sorted((item for item in items if all(
    item["source"].casefold() == term if term in ("custom", "coot") else
    term in (item["path"] + " " + item["label"]).casefold() for term in terms)),
    key=lambda item: (not item["label"].casefold().startswith(phrase), item["source"] != "Custom",
                      item["path"].casefold(), item["label"].casefold()))


def _native_menu_action_bridge():
  bridge = globals().get("_NATIVE_MENU_ACTION_BRIDGE")
  if bridge is None:
    import ctypes
    # Use the already loaded libraries, never a second GTK from another prefix.
    # Unlike GI/CDLL, PyDLL retains the GIL while Coot's C++ action re-enters Python.
    bridge = ctypes.PyDLL(None)
    pointer = ctypes.c_void_p
    signatures = {
      "PyCapsule_GetPointer": ([ctypes.py_object, ctypes.c_char_p], pointer),
      "gtk_widget_activate_action_variant": ([pointer, ctypes.c_char_p, pointer], ctypes.c_int),
      "g_variant_type_new": ([ctypes.c_char_p], pointer),
      "g_variant_type_free": ([pointer], None),
      "g_variant_parse": ([pointer, ctypes.c_char_p, pointer, pointer, pointer], pointer),
      "g_variant_take_ref": ([pointer], pointer),
      "g_variant_unref": ([pointer], None),
    }
    for name, (arguments, result) in signatures.items():
      function = getattr(bridge, name)
      function.argtypes, function.restype = arguments, result
    globals()["_NATIVE_MENU_ACTION_BRIDGE"] = bridge
  return bridge


def _activate_native_menu_action(widget, name, target=None):
  bridge = _native_menu_action_bridge()
  widget_pointer = bridge.PyCapsule_GetPointer(widget.__gpointer__, None)
  parameter = None
  if target is not None:
    # Reconstruct through the public GVariant API, not a private boxed-object layout.
    variant_type = bridge.g_variant_type_new(target.get_type_string().encode("ascii"))
    try:
      parameter = bridge.g_variant_parse(variant_type, target.print_(True).encode("utf-8"), None, None, None)
    finally:
      bridge.g_variant_type_free(variant_type)
    if not parameter:
      raise ValueError("The native menu action's parameter could not be converted safely.")
    bridge.g_variant_take_ref(parameter)
  try:
    return bool(bridge.gtk_widget_activate_action_variant(widget_pointer, name.encode("utf-8"), parameter))
  finally:
    if parameter:
      bridge.g_variant_unref(parameter)


def _install_custom_menu_search(menu):
  """Search Custom and Coot actions inside the existing Custom popover."""
  gui = _coot_gui_module()
  api = getattr(gui, "coot_gui_api", None)
  glib = _coot_gui_repository_module("GLib")
  if api is None or glib is None:
    return None
  toolbar = api.main_toolbar()
  button = toolbar.get_first_child()
  while button is not None:
    if isinstance(button, Gtk.MenuButton):
      popover = button.get_popover()
      if isinstance(popover, Gtk.PopoverMenu) and popover.get_menu_model() == menu:
        break
    button = button.get_next_sibling()
  if button is None:
    return None
  original_popover = popover
  app = api.application()
  body, root = Gio.Menu.new(), Gio.Menu.new()
  slot = Gio.MenuItem.new(None, None)
  slot.set_attribute_value("custom", glib.Variant("s", "coot-trimmings-search"))
  root.append_item(slot)
  root.append_section(None, body)
  # Cascading native submenus avoid GTK's duplicate sliding-page names when
  # the normal menu is removed/restored as a search begins or is cleared.
  popover = Gtk.PopoverMenu.new_from_model_full(root, Gtk.PopoverMenuFlags.NESTED)
  header = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6)
  for side in ("start", "end", "top", "bottom"):
    getattr(header, "set_margin_" + side)(6)
  search = Gtk.SearchEntry()
  search.set_placeholder_text("Search Custom and Coot menus...")
  search.set_size_request(340, -1)
  if hasattr(search, "set_search_delay"):
    search.set_search_delay(0)
  header.append(search)
  results = Gtk.ListBox()
  results.set_selection_mode(Gtk.SelectionMode.SINGLE)
  results.set_activate_on_single_click(True)
  scroller = Gtk.ScrolledWindow()
  scroller.set_policy(Gtk.PolicyType.NEVER, Gtk.PolicyType.AUTOMATIC)
  scroller.set_propagate_natural_height(True)
  scroller.set_max_content_height(360)
  scroller.set_child(results)
  header.append(scroller)
  empty = Gtk.Label(label="No matching menu items.")
  empty.add_css_class("dim-label")
  empty.set_margin_top(8)
  empty.set_margin_bottom(8)
  header.append(empty)
  state = {"search": search, "results": results, "body": body, "root": root, "popover": popover,
           "rows": [], "launching": False, "filtering": None}
  def live_items():
    get_native = getattr(api, "main_menumodel", None)
    native = get_native() if callable(get_native) else None
    return _menu_search_items(menu, native)
  def enabled(item):
    prefix, _dot, name = item["action"].partition(".")
    group = app if prefix == "app" else toolbar.get_root() if prefix == "win" else None
    lookup = getattr(group, "lookup_action", None)
    action = lookup(name) if callable(lookup) else None
    return bool(action and action.get_enabled())
  def refresh(*_args):
    child = results.get_first_child()
    while child is not None:
      results.remove(child)
      child = results.get_first_child()
    state["rows"] = []
    query = search.get_text().strip()
    if bool(query) != state["filtering"]:
      body.remove_all()
      if not query:
        body.append_section(None, menu)
      state["filtering"] = bool(query)
    if query:
      for item in _menu_search_matches(live_items(), query):
        row = Gtk.ListBoxRow()
        row.set_sensitive(enabled(item))
        labels = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=2)
        labels.set_margin_top(5)
        labels.set_margin_bottom(5)
        for text, dim in ((item["label"], False), (item["path"] or item["source"], True)):
          label = Gtk.Label(label=text)
          label.set_xalign(0)
          label.set_wrap(True)
          label.set_max_width_chars(52)
          if dim:
            label.add_css_class("dim-label")
          labels.append(label)
        labels.set_hexpand(True)
        content = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=12)
        content.append(labels)
        badge = Gtk.Label(label=item["source"])
        badge.add_css_class("dim-label")
        badge.set_valign(Gtk.Align.START)
        badge.set_margin_top(5)
        badge.set_tooltip_text("Coot-native menu action" if item["source"] == "Coot" else "coot_trimmings command")
        content.append(badge)
        row.set_child(content)
        row.set_tooltip_text(" > ".join(part for part in (item["source"], item["path"], item["label"]) if part))
        results.append(row)
        state["rows"].append((row, item))
      first = next((row for row, item in state["rows"] if row.get_sensitive()), None)
      results.select_row(first)
    scroller.set_visible(bool(query and state["rows"]))
    empty.set_visible(bool(query and not state["rows"]))
  def activate(_list=None, row=None):
    if state["launching"]:
      return
    row = row or results.get_selected_row()
    item = next((item for candidate, item in state["rows"] if candidate == row), None)
    if item is None or not enabled(item):
      return
    # Recent-file action names are reused. Never dispatch a result whose label
    # or target no longer matches the live menu entry behind that action name.
    if item not in live_items():
      refresh()
      add_status_bar_text("Menu changed; please select the search result again.")
      return
    state["launching"] = True
    popover.popdown()
    def run():
      state["launching"] = False
      try:
        if not enabled(item) or item not in live_items():
          add_status_bar_text("Menu changed; please select the search result again.")
          return False
        activated = (_activate_native_menu_action(button, item["action"], item["target"])
                     if item["source"] == "Coot" else button.activate_action(item["action"], item["target"]))
        if not activated:
          info_dialog("This menu action is no longer available.")
      except Exception as error:
        traceback.print_exc()
        info_dialog(f"Menu action failed: {error}")
      return False
    glib.idle_add(run)
  def key_pressed(_controller, keyval, _keycode, _modifiers):
    if keyval in (65362, 65364) and search.get_text().strip():
      rows = [row for row, item in state["rows"] if row.get_sensitive()]
      if rows:
        selected = results.get_selected_row()
        index = rows.index(selected) if selected in rows else (-1 if keyval == 65364 else 0)
        row = rows[(index + (1 if keyval == 65364 else -1)) % len(rows)]
        results.select_row(row)
        # Focusing the row scrolls it into view; keep typing in the search box.
        row.grab_focus()
        search.grab_focus()
      return True
    if keyval in (65293, 65421) and search.get_text().strip():
      activate()
      return True
    if keyval == 65307:
      if search.get_text():
        search.set_text("")
        refresh()
      else:
        popover.popdown()
      return True
    return False
  def opened(*_args):
    search.set_text("")
    refresh()
    def focus():
      if popover.get_mapped():
        search.grab_focus()
      return False
    glib.idle_add(focus)
  controller = Gtk.EventControllerKey()
  controller.set_propagation_phase(Gtk.PropagationPhase.CAPTURE)
  controller.connect("key-pressed", key_pressed)
  search.add_controller(controller)
  search.connect("search-changed", refresh)
  search.connect("activate", lambda *_: activate())
  results.connect("row-activated", activate)
  popover.connect("map", opened)
  state.update(refresh=refresh, activate=activate, key_pressed=key_pressed)
  try:
    if not popover.add_child(header, "coot-trimmings-search"):
      raise RuntimeError("Coot could not embed the Custom menu search field.")
    refresh()
    button.set_popover(popover)
  except Exception:
    button.set_popover(original_popover)
    raise
  button.custom_menu_search = state
  return state


def _toolbar_parse_call(text):
  """Accept one public function call, never statements or unrestricted eval."""
  if not isinstance(text, str) or not text.strip() or len(text) > 4000:
    raise ValueError("Enter one function call, for example cylinder_refine().")
  try:
    tree = ast.parse(text.strip(), mode="eval")
  except (SyntaxError, ValueError) as error:
    raise ValueError("Enter a Python function call, not a statement or code block.") from error
  allowed = (ast.Expression, ast.Call, ast.Name, ast.Load, ast.Attribute, ast.Constant,
             ast.List, ast.Tuple, ast.Dict, ast.UnaryOp, ast.UAdd, ast.USub,
             ast.Subscript, ast.Slice, ast.keyword)
  nodes = list(ast.walk(tree))
  if len(nodes) > 150 or not isinstance(tree.body, ast.Call):
    raise ValueError("Enter a single, reasonably short function call.")
  for node in nodes:
    if not isinstance(node, allowed):
      raise ValueError("Use calls, literal arguments, keywords and indexing; no operators or code blocks.")
    if isinstance(node, ast.Name) and node.id.startswith("_"):
      raise ValueError("Private names are not supported in toolbar calls.")
    if isinstance(node, ast.Attribute) and (not isinstance(node.value, ast.Name)
        or node.value.id not in ("coot", "coot_utils", "coot_fitting", "coot_gui") or node.attr.startswith("_")):
      raise ValueError("Qualified functions must be public members of coot, coot_utils, coot_fitting or coot_gui.")
    if isinstance(node, ast.Call) and not isinstance(node.func, (ast.Name, ast.Attribute)):
      raise ValueError("Call a named function, not a computed callable.")
    if isinstance(node, ast.Call) and len({kw.arg for kw in node.keywords}) != len(node.keywords):
      raise ValueError("Keyword arguments must not be repeated.")
    if isinstance(node, ast.keyword) and (node.arg is None or node.arg.startswith("_")):
      raise ValueError("Expanded or private keyword arguments are not supported.")
    if isinstance(node, ast.Dict) and any(key is None for key in node.keys):
      raise ValueError("Expanded dictionaries are not supported.")
  return tree.body


def _toolbar_eval_call(node, execute=False):
  """Resolve names while editing/loading; execute calls only after a click."""
  if isinstance(node, ast.Constant):
    return node.value
  if isinstance(node, ast.Name):
    if node.id not in globals():
      raise ValueError(f"Function or value unavailable in this session: {node.id}")
    return globals()[node.id]
  if isinstance(node, ast.Attribute):
    module = globals().get(node.value.id)
    if module is None and node.value.id == "coot_gui":
      module = _coot_gui_module()
    if module is None or not hasattr(module, node.attr):
      raise ValueError(f"API unavailable in this build: {node.value.id}.{node.attr}")
    return getattr(module, node.attr)
  if isinstance(node, ast.Call):
    function = _toolbar_eval_call(node.func, execute)
    if not callable(function):
      raise ValueError("The selected name is not a function.")
    args = [_toolbar_eval_call(arg, execute) for arg in node.args]
    kwargs = {kw.arg: _toolbar_eval_call(kw.value, execute) for kw in node.keywords}
    return function(*args, **kwargs) if execute else None
  if isinstance(node, (ast.List, ast.Tuple)):
    values = [_toolbar_eval_call(value, execute) for value in node.elts]
    return tuple(values) if isinstance(node, ast.Tuple) else values
  if isinstance(node, ast.Dict):
    keys = [_toolbar_eval_call(key, execute) for key in node.keys]
    values = [_toolbar_eval_call(value, execute) for value in node.values]
    return dict(zip(keys, values)) if execute else None
  if isinstance(node, ast.UnaryOp):
    value = _toolbar_eval_call(node.operand, execute)
    return (-value if isinstance(node.op, ast.USub) else +value) if execute else None
  if isinstance(node, ast.Slice):
    values = [_toolbar_eval_call(value, execute) if value is not None else None
              for value in (node.lower, node.upper, node.step)]
    return slice(*values) if execute else None
  if isinstance(node, ast.Subscript):
    value, index = _toolbar_eval_call(node.value, execute), _toolbar_eval_call(node.slice, execute)
    return value[index] if execute else None
  raise ValueError("Unsupported toolbar call expression.")


def _toolbar_normalize_buttons(buttons):
  if not isinstance(buttons, list) or len(buttons) > 100:
    raise ValueError("Toolbar preferences must contain a list of at most 100 buttons.")
  result, ids = [], set()
  for button in buttons:
    if not isinstance(button, dict):
      raise ValueError("Invalid toolbar button definition.")
    row = {}
    for field, limit in (("id", 80), ("label", 100), ("tooltip", 500), ("icon", 250)):
      value = button.get(field, "")
      if not isinstance(value, str) or len(value) > limit:
        raise ValueError(f"Invalid toolbar {field}.")
      row[field] = value.strip()
    if not row["id"] or row["id"] in ids or not row["label"]:
      raise ValueError("Each button needs a label and a unique identifier.")
    ids.add(row["id"])
    for field in ("visible", "show_label"):
      row[field] = button.get(field, True)
      if type(row[field]) is not bool:
        raise ValueError(f"Invalid toolbar {field} setting.")
    icon = row["icon"]
    if icon and (not icon.startswith(("theme:", "coot:", "custom:")) or not icon.split(":", 1)[1]
                 or any(char in icon for char in ("/", "\\", "\0")) or ".." in icon):
      raise ValueError("Choose or import an icon using the icon picker.")
    row["kind"] = button.get("kind")
    if row["kind"] == "call":
      row["call"] = button.get("call", "")
      _toolbar_parse_call(row["call"])
    elif row["kind"] == "keybinding":
      action = button.get("action")
      if not isinstance(action, str) or not action.strip() or len(action) > 500:
        raise ValueError("Choose a keybinding action.")
      row["action"] = action
    elif row["kind"] == "menu":
      source = button.get("source")
      if not isinstance(source, dict) or not isinstance(source.get("path"), str) or not isinstance(source.get("label"), str):
        raise ValueError("Choose an item from the Custom menu.")
      if not source["label"] or len(source["path"]) > 1000 or len(source["label"]) > 500:
        raise ValueError("Invalid menu item identity.")
      target = source.get("target")
      if target is not None and (not isinstance(target, dict) or set(target) != {"type", "value"}
                                 or not isinstance(target["type"], str) or len(target["type"]) > 100):
        raise ValueError("Invalid menu action target.")
      row["source"] = {"path": source["path"], "label": source["label"], "target": target}
    else:
      raise ValueError("Choose a menu action, keybinding action or function call.")
    result.append(row)
  # Check that preferences really are JSON data, never widgets or callables.
  json.dumps(result, allow_nan=False)
  return result


def _toolbar_state():
  state = _CUSTOM_TOOLBAR_STATE
  if "buttons" not in state:
    state.update(buttons=[], menu=None, group=None, window=None, load_error=None)
    path = os.path.join(_coot_trimmings_dir(), CUSTOM_TOOLBAR_FILENAME)
    try:
      with open(path, encoding="utf-8") as handle:
        source = handle.read(1024*1024+1)
      if len(source) > 1024*1024:
        raise ValueError("Toolbar preferences file is too large.")
      data = json.loads(source)
      if not isinstance(data, dict) or data.get("version") != 1:
        raise ValueError("Unsupported toolbar preferences version.")
      state["buttons"] = _toolbar_normalize_buttons(data.get("buttons"))
    except FileNotFoundError:
      pass
    except Exception as error:
      state["load_error"] = (f"Could not read {path}: {error}. File left unchanged; "
                             "repair or remove it and restart Coot before saving.")
      print(state["load_error"])
  return state


def _toolbar_commit(buttons):
  state = _toolbar_state()
  if state["load_error"]:
    raise ValueError(state["load_error"])
  buttons = _toolbar_normalize_buttons(buttons)
  path = os.path.join(_coot_trimmings_dir(), CUSTOM_TOOLBAR_FILENAME)
  temporary = None
  try:
    fd, temporary = tempfile.mkstemp(prefix=".coot_toolbar_", suffix=".json", dir=os.path.dirname(path))
    with os.fdopen(fd, "w", encoding="utf-8") as handle:
      json.dump({"version": 1, "buttons": buttons}, handle, indent=2, ensure_ascii=True, allow_nan=False)
      handle.write("\n")
    os.replace(temporary, path)
    temporary = None
  finally:
    if temporary is not None:
      os.unlink(temporary)
  state["buttons"] = buttons
  _toolbar_render()


def _toolbar_menu_source(item):
  target = item["target"]
  return {"path": item["path"], "label": item["label"],
          "target": None if target is None else {"type": target.get_type_string(),
              "value": json.loads(json.dumps(target.unpack(), allow_nan=False))}}


def _toolbar_menu_items(query=""):
  menu = _toolbar_state()["menu"]
  if menu is None:
    return []
  items = _custom_menu_search_matches(menu, query) if query.strip() else _custom_menu_search_items(menu)
  # A recent-file slot can point at a different file next session. Pin commands,
  # not history slots; favourite monomers are resolved by their current label.
  return [item for item in items if not item["path"].startswith(("Recent maps", "Recent models"))]


def _toolbar_resolve_menu(source):
  matches = [item for item in _toolbar_menu_items() if _toolbar_menu_source(item) == source]
  if len(matches) != 1:
    raise ValueError("Menu item unavailable or ambiguous: " + source["path"] + " > " + source["label"])
  return matches[0]


def _toolbar_keybinding_items(query=""):
  terms = query.casefold().split()
  items = list(globals().get("REGISTERED_KEYBINDING_ACTIONS", {}).values())
  return sorted((item for item in items if all(term in
      (item["name"] + " " + _custom_key_binding_label(item["key"])).casefold() for term in terms)),
      key=lambda item: item["name"].casefold())


def _toolbar_resolve_keybinding(name):
  matches = [item for item in _toolbar_keybinding_items() if item["name"] == name]
  if len(matches) != 1:
    raise ValueError("Keybinding action unavailable or ambiguous: " + name)
  return matches[0]["callback"]


def _toolbar_problem(button):
  try:
    if button["kind"] == "call":
      _toolbar_eval_call(_toolbar_parse_call(button["call"]), execute=False)
    elif button["kind"] == "keybinding":
      _toolbar_resolve_keybinding(button["action"])
    else:
      item = _toolbar_resolve_menu(button["source"])
      action = _coot_gui_module().coot_gui_api.application().lookup_action(item["action"].removeprefix("app."))
      if action is None or not action.get_enabled():
        raise ValueError("Menu action is not currently available.")
  except Exception as error:
    return str(error)
  return ""


def _toolbar_icon_theme(widget):
  return Gtk.IconTheme.get_for_display(widget.get_display())


def _toolbar_custom_icon_directory():
  return os.path.join(_coot_trimmings_dir(), "coot_trimmings_toolbar_icons")


def _toolbar_custom_icon_name(token):
  if not isinstance(token, str) or not token.startswith("custom:"):
    raise ValueError("Only imported custom icons can be renamed or removed.")
  name = token.split(":", 1)[1]
  stem, suffix = os.path.splitext(name)
  if (len(stem) != 64 or any(char not in "0123456789abcdef" for char in stem)
      or suffix not in (".png", ".svg", ".jpg", ".jpeg", ".xpm", ".webp")):
    raise ValueError("Invalid custom icon filename.")
  return name


def _toolbar_validate_icon_names(names):
  if not isinstance(names, dict):
    raise ValueError("Invalid custom icon names.")
  for filename, label in names.items():
    if not isinstance(filename, str):
      raise ValueError("Invalid custom icon filename.")
    _toolbar_custom_icon_name("custom:" + filename)
    if not isinstance(label, str) or not label.strip() or len(label) > 100 or not label.isprintable():
      raise ValueError("Icon names must contain 1-100 printable characters.")
  return names


def _toolbar_icon_names():
  path = os.path.join(_toolbar_custom_icon_directory(), ".names.json")
  try:
    with open(path, encoding="utf-8") as handle:
      source = handle.read(256*1024+1)
  except FileNotFoundError:
    return {}
  if len(source) > 256*1024:
    raise ValueError("Custom icon names file is too large.")
  data = json.loads(source)
  if not isinstance(data, dict) or data.get("version") != 1:
    raise ValueError("Unsupported custom icon names version.")
  return _toolbar_validate_icon_names(data.get("names"))


def _toolbar_write_icon_names(names):
  names = _toolbar_validate_icon_names(names)
  data = json.dumps({"version": 1, "names": names}, indent=2, ensure_ascii=True) + "\n"
  if len(data) > 256*1024:
    raise ValueError("Custom icon names file is too large.")
  directory = _toolbar_custom_icon_directory()
  os.makedirs(directory, exist_ok=True)
  temporary = None
  try:
    fd, temporary = tempfile.mkstemp(prefix=".names_", suffix=".json", dir=directory)
    with os.fdopen(fd, "w", encoding="utf-8") as handle:
      handle.write(data)
    os.replace(temporary, os.path.join(directory, ".names.json"))
    temporary = None
  finally:
    if temporary is not None:
      os.unlink(temporary)


def _toolbar_rename_icon(token, label):
  filename = _toolbar_custom_icon_name(token)
  if not isinstance(label, str):
    raise ValueError("Enter an icon name.")
  label = label.strip()
  _toolbar_validate_icon_names({filename: label})
  if not os.path.isfile(os.path.join(_toolbar_custom_icon_directory(), filename)):
    raise ValueError("This custom icon is no longer available.")
  names = _toolbar_icon_names()
  names[filename] = label
  _toolbar_write_icon_names(names)
  return label


def _toolbar_import_icon(path):
  suffix = os.path.splitext(path)[1].lower()
  if suffix not in (".png", ".svg", ".jpg", ".jpeg", ".xpm", ".webp"):
    raise ValueError("Choose a PNG, SVG, JPEG, XPM or WebP image.")
  with open(path, "rb") as handle:
    data = handle.read(5*1024*1024+1)
  if not data or len(data) > 5*1024*1024:
    raise ValueError("Icon images must be nonempty and smaller than 5 MB.")
  pixbuf = _coot_gui_repository_module("GdkPixbuf")
  if pixbuf is None:
    raise ValueError("Image loading is unavailable in this Coot build.")
  loader = pixbuf.PixbufLoader.new()
  loader.set_size(48, 48)
  try:
    loader.write(data)
    loader.close()
    if loader.get_pixbuf() is None:
      raise ValueError("This image could not be loaded.")
  except Exception as error:
    raise ValueError("This image could not be loaded: " + str(error)) from error
  name = hashlib.sha256(data).hexdigest() + suffix
  names = _toolbar_icon_names()
  directory = _toolbar_custom_icon_directory()
  os.makedirs(directory, exist_ok=True)
  temporary = None
  try:
    fd, temporary = tempfile.mkstemp(prefix=".icon_", dir=directory)
    with os.fdopen(fd, "wb") as handle:
      handle.write(data)
    os.replace(temporary, os.path.join(directory, name))
    temporary = None
  finally:
    if temporary is not None:
      os.unlink(temporary)
  if name not in names:
    names[name] = "".join(char if char.isprintable() else " " for char in os.path.basename(path)).strip()[:100] or "Imported icon"
    _toolbar_write_icon_names(names)
  return "custom:" + name


def _toolbar_remove_custom_icon(token):
  name = _toolbar_custom_icon_name(token)
  try:
    os.unlink(os.path.join(_toolbar_custom_icon_directory(), name))
  except FileNotFoundError:
    pass
  # Keep its label so re-importing the same image restores the searchable name.
  _toolbar_render()


def _toolbar_add_icon_drop(widget, import_path, show_error):
  gdk = _coot_gui_repository_module("Gdk")
  if gdk is None or not hasattr(gdk, "FileList") or not hasattr(Gtk, "DropTarget"):
    return None
  target = Gtk.DropTarget.new(gdk.FileList.__gtype__, gdk.DragAction.COPY)
  def dropped(_target, value, _x, _y):
    try:
      files = value.get_files()
      if len(files) != 1:
        raise ValueError("Drop one image file at a time.")
      path = files[0].get_path()
      if not path:
        raise ValueError("Drop an image saved on this computer.")
      return bool(import_path(path))
    except Exception as error:
      show_error(str(error))
      return False
  target.connect("drop", dropped)
  widget.add_controller(target)
  return target


def _toolbar_add_file_preview(chooser):
  glib = _coot_gui_repository_module("GLib")
  gio = _coot_gui_repository_module("Gio")
  pixbuf = _coot_gui_repository_module("GdkPixbuf")
  if (glib is None or gio is None or pixbuf is None
      or not hasattr(pixbuf.Pixbuf, "new_from_stream_at_scale_async")):
    return None
  pane = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=8)
  pane.set_vexpand(True)
  for side in ("start", "end", "top", "bottom"):
    getattr(pane, "set_margin_"+side)(10)
  details = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=12)
  details.set_hexpand(True)
  details.append(Gtk.Label(label="Selected image", xalign=0))
  description = Gtk.Label(label="Select an image to preview.", xalign=0)
  description.set_wrap(True)
  description.set_max_width_chars(45)
  description.set_hexpand(True)
  details.append(description)
  pane.append(details)
  # Gtk.Image applies icon sizing. Picture scales the cached preview to its
  # allocation, including on resize, without decoding the selected file again.
  image = Gtk.Picture()
  image.set_can_shrink(True)
  if hasattr(image, "set_content_fit") and hasattr(Gtk, "ContentFit"):
    image.set_content_fit(Gtk.ContentFit.CONTAIN)
  else:
    image.set_keep_aspect_ratio(True)
  image.set_hexpand(True)
  image.set_vexpand(True)
  image.set_size_request(-1, 180)
  pane.append(image)
  chooser.get_content_area().append(pane)
  state = {"closed": False, "source": None, "path": None, "started": False,
           "cancellable": None, "loads": 0}
  def cancel_load():
    if state["cancellable"] is not None:
      state["cancellable"].cancel()
      state["cancellable"] = None
  def stop(*_args):
    state["closed"] = True
    if state["source"] is not None:
      glib.source_remove(state["source"])
      state["source"] = None
    cancel_load()
    return False
  def close_stream(stream):
    if stream is not None:
      try:
        stream.close(None)
      except Exception:
        pass
  def start_load(file, path):
    request = gio.Cancellable.new()
    state["cancellable"] = request
    def current():
      return not state["closed"] and state["cancellable"] is request and not request.is_cancelled()
    def failed(error):
      if current():
        description.set_text(os.path.basename(path) + "\nPreview unavailable: " + str(error))
    def opened(source, result, *_args):
      stream = None
      try:
        stream = source.read_finish(result)
        if not current():
          close_stream(stream)
          return
        def decoded(_source, result, *_args):
          try:
            thumbnail = pixbuf.Pixbuf.new_from_stream_finish(result)
            if current():
              image.set_pixbuf(thumbnail)
              description.set_text(os.path.basename(path))
          except Exception as error:
            failed(error)
          finally:
            close_stream(stream)
        state["loads"] += 1
        pixbuf.Pixbuf.new_from_stream_at_scale_async(stream, 512, 384, True, request, decoded)
      except Exception as error:
        close_stream(stream)
        failed(error)
    def inspected(source, result, *_args):
      try:
        info = source.query_info_finish(result)
        if not current():
          return
        if info.get_file_type() != gio.FileType.REGULAR:
          description.set_text("Select an image file to preview.")
        elif not 0 < info.get_size() <= 5*1024*1024:
          description.set_text("Preview skipped: icon images must be nonempty and at most 5 MB.")
        else:
          file.read_async(glib.PRIORITY_DEFAULT, request, opened)
      except Exception as error:
        failed(error)
    file.query_info_async("standard::size,standard::type", gio.FileQueryInfoFlags.NONE,
                          glib.PRIORITY_DEFAULT, request, inspected)
  def tick():
    if state["closed"]:
      state["source"] = None
      return False
    if not chooser.get_mapped():
      return True
    try:
      file = chooser.get_file()
      path = file.get_path() if file is not None else None
      if path != state["path"]:
        cancel_load()
        state.update(path=path, started=False)
        image.set_paintable(None)
        description.set_text(os.path.basename(path) if path else "Select an image to preview.")
        return True
      if path and not state["started"]:
        state["started"] = True
        if os.path.splitext(path)[1].lower() in (".png", ".svg", ".jpg", ".jpeg", ".xpm", ".webp"):
          start_load(file, path)
        else:
          description.set_text("Select a supported image file to preview.")
    except Exception:
      description.set_text("Preview unavailable.")
      stop()
      return False
    return True
  # GTK4 has no public selection-changed signal on FileChooser. Only poll while
  # this browser is open; two equal samples debounce navigation without file I/O.
  state["source"] = glib.timeout_add(250, tick)
  chooser.connect("unrealize", stop)
  chooser.toolbar_preview_widgets = {"image": image, "description": description,
                                     "state": state, "stop": stop, "tick": tick}
  return stop


def _toolbar_icon_image(token, widget):
  if not token:
    return None
  theme = _toolbar_icon_theme(widget)
  kind, name = token.split(":", 1)
  image = None
  if kind == "theme" and theme.has_icon(name):
    image = Gtk.Image.new_from_icon_name(name)
  elif kind == "coot":
    for directory in theme.get_search_path():
      path = os.path.join(directory, name)
      if os.path.isfile(path):
        image = Gtk.Image.new_from_file(path)
        break
  elif kind == "custom":
    path = os.path.join(_toolbar_custom_icon_directory(), name)
    if os.path.isfile(path):
      image = Gtk.Image.new_from_file(path)
  if image is not None:
    image.set_pixel_size(24)
  return image


def _toolbar_icon_catalog(widget):
  theme = _toolbar_icon_theme(widget)
  result = {"theme:"+name: name for name in theme.get_icon_names()}
  directory = _toolbar_custom_icon_directory()
  try:
    names = _toolbar_icon_names()
  except Exception as error:
    names = {}
    print("Could not read custom icon names; file left unchanged:", error)
  if os.path.isdir(directory):
    for name in sorted(os.listdir(directory)):
      if not name.startswith(".") and os.path.isfile(os.path.join(directory, name)):
        result["custom:"+name] = "Custom: " + names.get(name, name[:12] + os.path.splitext(name)[1])
  for directory in theme.get_search_path():
    try:
      for name in os.listdir(directory):
        if name.lower().endswith((".svg", ".png", ".xpm")) and os.path.isfile(os.path.join(directory, name)):
          result.setdefault("coot:"+name, "Coot: " + name)
    except OSError:
      continue
  return sorted(result.items(), key=lambda item: (0 if item[0].startswith("custom:") else
      1 if item[0].startswith("coot:") else 2, item[1].casefold()))


def _toolbar_activate(button_id, widget):
  button = next((row for row in _toolbar_state()["buttons"] if row["id"] == button_id), None)
  if button is None:
    return
  try:
    problem = _toolbar_problem(button)
    if problem:
      raise ValueError(problem)
    if button["kind"] == "call":
      return _toolbar_eval_call(_toolbar_parse_call(button["call"]), execute=True)
    if button["kind"] == "keybinding":
      return _toolbar_resolve_keybinding(button["action"])()
    item = _toolbar_resolve_menu(button["source"])
    if not widget.activate_action(item["action"], item["target"]):
      raise ValueError("This menu action is no longer available.")
  except Exception as error:
    traceback.print_exc()
    info_dialog(f"Toolbar button '{button['label']}' could not run:\n\n{error}")


def _toolbar_render():
  state = _toolbar_state()
  group = state.get("group")
  if group is None:
    return
  child = group.get_first_child()
  while child is not None:
    group.remove(child)
    child = group.get_first_child()
  for definition in state["buttons"]:
    if not definition["visible"]:
      continue
    button = Gtk.Button()
    content = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=4)
    image = None
    try:
      image = _toolbar_icon_image(definition["icon"], group)
    except Exception as error:
      print("Toolbar icon unavailable:", error)
    if image is not None:
      content.append(image)
    if definition["show_label"] or image is None:
      content.append(Gtk.Label(label=definition["label"]))
    button.set_child(content)
    problem = _toolbar_problem(definition)
    tooltip = definition["tooltip"] or definition["label"]
    if definition["icon"] and image is None:
      tooltip += "\nIcon unavailable; showing text."
    button.set_tooltip_text(tooltip + ("\nUnavailable: " + problem if problem else ""))
    button.set_sensitive(not problem)
    button.custom_toolbar_id = definition["id"]
    button.connect("clicked", lambda widget, button_id=definition["id"]: _toolbar_activate(button_id, widget))
    group.append(button)
  group.set_visible(any(button["visible"] for button in state["buttons"]))


def _toolbar_choose_icon(parent, selected, initial="", library_changed=None):
  window = Gtk.Window(title="Choose toolbar icon")
  window.set_transient_for(parent)
  window.set_destroy_with_parent(True)
  window.set_default_size(540, 440)
  box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=8)
  for side in ("start", "end", "top", "bottom"):
    getattr(box, "set_margin_"+side)(10)
  search = Gtk.SearchEntry()
  search.set_placeholder_text("Search icons...")
  box.append(search)
  scroller = Gtk.ScrolledWindow()
  scroller.set_vexpand(True)
  grid = Gtk.FlowBox()
  grid.set_selection_mode(Gtk.SelectionMode.NONE)
  grid.set_max_children_per_line(6)
  scroller.set_child(grid)
  box.append(scroller)
  preview_row = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
  preview = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
  preview.set_hexpand(True)
  use_icon = Gtk.Button(label="Use icon")
  use_icon.set_sensitive(False)
  preview_row.append(preview)
  preview_row.append(use_icon)
  box.append(preview_row)
  name_row = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
  name_row.append(Gtk.Label(label="Icon name"))
  icon_name = Gtk.Entry()
  icon_name.set_max_length(100)
  icon_name.set_hexpand(True)
  rename_icon = Gtk.Button(label="Rename")
  name_row.append(icon_name)
  name_row.append(rename_icon)
  box.append(name_row)
  library_row = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
  import_icon = Gtk.Button(label="Import image...")
  import_icon.set_tooltip_text("Copy your own image into the saved toolbar icons folder")
  remove_icon = Gtk.Button(label="Remove custom icon")
  remove_icon.set_sensitive(False)
  library_row.append(import_icon)
  library_row.append(remove_icon)
  box.append(library_row)
  status = Gtk.Label(xalign=0)
  status.set_wrap(True)
  box.append(status)
  footer = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
  count = Gtk.Label()
  more = Gtk.Button(label="Load more")
  no_icon = Gtk.Button(label="Text only")
  close = Gtk.Button(label="Cancel")
  for widget in (count, more, no_icon, close):
    footer.append(widget)
  box.append(footer)
  window.set_child(box)
  catalog = _toolbar_icon_catalog(parent)
  state = {"limit": 96, "chooser": None, "closed": False, "parent_handler": None,
           "selected": "", "confirm_remove": False, "tiles": []}
  def dismiss_chooser():
    chooser = state["chooser"]
    state["chooser"] = None
    if chooser is not None:
      stop_preview = getattr(chooser, "toolbar_stop_preview", None)
      if stop_preview is not None:
        stop_preview()
      chooser.destroy()
  def cleanup(*_args):
    state["closed"] = True
    dismiss_chooser()
    if state["parent_handler"] is not None:
      parent.disconnect(state["parent_handler"])
      state["parent_handler"] = None
    return False
  def close_picker(*_args):
    cleanup()
    window.destroy()
  def choose(token):
    if state["closed"]:
      return
    selected(token)
    close_picker()
  def preview_icon(token):
    if state["closed"]:
      return
    state.update(selected=token, confirm_remove=False)
    child = preview.get_first_child()
    while child is not None:
      preview.remove(child)
      child = preview.get_first_child()
    image = _toolbar_icon_image(token, parent)
    if image is not None:
      image.set_pixel_size(48)
      preview.append(image)
    name = dict(catalog).get(token, "Choose an icon to preview")
    preview.append(Gtk.Label(label=name if len(name) <= 35 else name[:32]+"...", xalign=0))
    preview.set_tooltip_text(name)
    use_icon.set_sensitive(bool(token))
    remove_icon.set_label("Remove custom icon")
    remove_icon.set_sensitive(token.startswith("custom:"))
    name_row.set_visible(token.startswith("custom:"))
    icon_name.set_text(name.removeprefix("Custom: ") if token.startswith("custom:") else "")
    for tile, candidate in state["tiles"]:
      if candidate == token:
        tile.add_css_class("suggested-action")
      else:
        tile.remove_css_class("suggested-action")
    status.set_text("")
  def render(reset=False):
    child = grid.get_first_child()
    while child is not None:
      grid.remove(child)
      child = grid.get_first_child()
    state["tiles"] = []
    terms = search.get_text().casefold().split()
    matching = [(token, label) for token, label in catalog if all(term in label.casefold() for term in terms)]
    for token, label in matching[:state["limit"]]:
      button = Gtk.Button()
      content = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=3)
      image = _toolbar_icon_image(token, parent)
      if image is not None:
        content.append(image)
      name = label.removeprefix("Coot: ").removeprefix("Custom: ")
      content.append(Gtk.Label(label=name if len(name) <= 13 else name[:10]+"..."))
      button.set_child(content)
      button.set_tooltip_text(label)
      if token == state["selected"]:
        button.add_css_class("suggested-action")
      button.connect("clicked", lambda _widget, token=token: preview_icon(token))
      grid.append(button)
      state["tiles"].append((button, token))
    count.set_text(f"{min(len(matching), state['limit'])} / {len(matching)} icons")
    more.set_visible(len(matching) > state["limit"])
  def load_more(*_args):
    state["limit"] += 96
    render()
  def import_path(path):
    if state["closed"]:
      return False
    try:
      token = _toolbar_import_icon(path)
      catalog[:] = _toolbar_icon_catalog(parent)
      search.set_text("")
      render()
      preview_icon(token)
      scroller.get_vadjustment().set_value(0)
      status.set_text("Added to custom icons. Click Use icon to apply it to the button.")
      return True
    except Exception as error:
      status.set_text(str(error))
      return False
  def remove_selected(*_args):
    token = state["selected"]
    if state["closed"] or not token.startswith("custom:"):
      return
    users = [button["label"] for button in _toolbar_state()["buttons"] if button["icon"] == token]
    if not state["confirm_remove"]:
      state["confirm_remove"] = True
      remove_icon.set_label("Confirm removal")
      warning = (f" Used by {len(users)} saved button(s); these will fall back to text." if users else "")
      status.set_text("Remove this image from custom icons?" + warning + " The original file is kept.")
      return
    try:
      _toolbar_remove_custom_icon(token)
      catalog[:] = _toolbar_icon_catalog(parent)
      preview_icon("")
      render()
      if library_changed is not None:
        library_changed()
      status.set_text("Removed from custom icons." + (" Buttons using it now show text." if users else ""))
    except Exception as error:
      status.set_text("Could not remove icon: " + str(error))
  def rename_selected(*_args):
    token = state["selected"]
    if state["closed"] or not token.startswith("custom:"):
      return
    try:
      label = _toolbar_rename_icon(token, icon_name.get_text())
      catalog[:] = _toolbar_icon_catalog(parent)
      search.set_text("")
      render()
      preview_icon(token)
      status.set_text("Icon renamed to '" + label + "'.")
    except Exception as error:
      status.set_text("Could not rename icon: " + str(error))
  def browse_image(*_args):
    if state["closed"]:
      return
    if state["chooser"] is not None:
      state["chooser"].present()
      return
    # Avoid the Quartz native chooser's asynchronous completion after its parent
    # disappears. A GTK dialog also respects the picker window's stacking order.
    chooser = Gtk.FileChooserDialog(title="Import toolbar icon", transient_for=window,
        modal=True, action=Gtk.FileChooserAction.OPEN)
    chooser.set_default_size(800, 640)
    chooser.set_destroy_with_parent(True)
    chooser.add_buttons("Cancel", Gtk.ResponseType.CANCEL, "Import", Gtk.ResponseType.ACCEPT)
    chooser.set_default_response(Gtk.ResponseType.ACCEPT)
    images = Gtk.FileFilter()
    images.set_name("Icon images")
    for suffix in ("png", "svg", "jpg", "jpeg", "xpm", "webp"):
      images.add_pattern("*."+suffix)
      images.add_pattern("*."+suffix.upper())
    chooser.add_filter(images)
    chooser.toolbar_stop_preview = _toolbar_add_file_preview(chooser)
    def response(dialog, response_id):
      if state["closed"] or state["chooser"] is not dialog:
        return
      chosen = dialog.get_file() if response_id == Gtk.ResponseType.ACCEPT else None
      path = chosen.get_path() if chosen is not None else None
      dismiss_chooser()
      window.present()
      if path:
        import_path(path)
    def chooser_closed(dialog):
      response(dialog, Gtk.ResponseType.CANCEL)
      return True
    chooser.connect("response", response)
    chooser.connect("close-request", chooser_closed)
    state["chooser"] = chooser
    chooser.present()
  window.connect("close-request", cleanup)
  window.connect("unrealize", cleanup)
  state["parent_handler"] = parent.connect("unrealize", close_picker)
  import_icon.connect("clicked", browse_image)
  search.connect("search-changed", lambda *_: render(True))
  more.connect("clicked", load_more)
  use_icon.connect("clicked", lambda *_: choose(state["selected"]) if state["selected"] else None)
  remove_icon.connect("clicked", remove_selected)
  rename_icon.connect("clicked", rename_selected)
  icon_name.connect("activate", rename_selected)
  no_icon.connect("clicked", lambda *_: choose(""))
  close.connect("clicked", close_picker)
  drop_target = _toolbar_add_icon_drop(box, import_path, status.set_text)
  if drop_target is not None:
    hint = Gtk.Label(label="Or drop an image file into this window.")
    hint.add_css_class("dim-label")
    box.append(hint)
  window.toolbar_icon_widgets = {"search": search, "grid": grid, "more": more,
      "text_only": no_icon, "render": render, "choose": choose, "state": state,
      "count": count, "import_icon": import_icon, "import_path": import_path, "status": status,
      "cancel": close, "drop_target": drop_target, "preview_icon": preview_icon,
      "use_icon": use_icon, "remove_icon": remove_icon, "preview": preview,
      "icon_name": icon_name, "rename_icon": rename_icon, "name_row": name_row}
  render()
  preview_icon(initial if initial in dict(catalog) else "")
  window.present()
  return window


def custom_toolbar_dialog():
  state = _toolbar_state()
  if state.get("window") is not None:
    state["window"].present()
    return state["window"]
  if globals().get("Gtk") is None or state.get("group") is None:
    info_dialog("Custom toolbar buttons require Coot's GTK toolbar.")
    return None
  _toolbar_render()
  window = Gtk.Window(title="Customize toolbar")
  window.set_default_size(760, 620)
  outer = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=10)
  for side in ("start", "end", "top", "bottom"):
    getattr(outer, "set_margin_"+side)(12)
  columns = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=14)
  columns.set_vexpand(True)
  outer.append(columns)
  left = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6)
  left.set_size_request(220, -1)
  left.append(Gtk.Label(label="Saved buttons", xalign=0))
  saved = Gtk.ListBox()
  saved.set_selection_mode(Gtk.SelectionMode.NONE)
  saved.set_activate_on_single_click(True)
  scroll = Gtk.ScrolledWindow()
  scroll.set_vexpand(True)
  scroll.set_child(saved)
  left.append(scroll)
  columns.append(left)
  right = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=8)
  right.set_hexpand(True)
  columns.append(right)
  controls = {}
  heading = Gtk.Label(label="Choose a button", xalign=0)
  heading.add_css_class("heading")
  right.append(heading)
  def field(label, widget):
    right.append(Gtk.Label(label=label, xalign=0))
    right.append(widget)
    return widget
  label = field("Button label", Gtk.Entry())
  tooltip = field("Tooltip (optional)", Gtk.Entry())
  mode = Gtk.ComboBoxText()
  for text in ("Custom menu action", "Function call", "Keybinding action"):
    mode.append_text(text)
  field("Action", mode)
  pages = Gtk.Stack()
  right.append(pages)
  menu_page = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=5)
  menu_search = Gtk.SearchEntry()
  menu_search.set_placeholder_text("Find a menu action...")
  menu_page.append(menu_search)
  choices = Gtk.ListBox()
  choices.set_selection_mode(Gtk.SelectionMode.SINGLE)
  menu_scroll = Gtk.ScrolledWindow()
  menu_scroll.set_policy(Gtk.PolicyType.NEVER, Gtk.PolicyType.AUTOMATIC)
  menu_scroll.set_size_request(-1, 170)
  menu_scroll.set_child(choices)
  menu_page.append(menu_scroll)
  pages.add_named(menu_page, "menu")
  call_page = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6)
  call = Gtk.Entry()
  call.add_css_class("monospace")
  call.set_placeholder_text("cylinder_refine()")
  call_page.append(call)
  hint = Gtk.Label(label="One API or script helper call, with arguments if needed.\n"
                         "Runs only when clicked and may change your model.", xalign=0)
  hint.set_wrap(True)
  hint.add_css_class("dim-label")
  call_page.append(hint)
  pages.add_named(call_page, "call")
  source_label = Gtk.Label(xalign=0)
  source_label.set_wrap(True)
  right.append(source_label)
  icon_row = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
  preview = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=4)
  choose_icon = Gtk.Button(label="Choose icon...")
  clear_icon = Gtk.Button(label="No icon")
  for widget in (preview, choose_icon, clear_icon):
    icon_row.append(widget)
  right.append(icon_row)
  show_label = Gtk.CheckButton(label="Show label beside icon")
  visible = Gtk.CheckButton(label="Show button in toolbar")
  right.append(show_label)
  right.append(visible)
  status = Gtk.Label(xalign=0)
  status.set_wrap(True)
  outer.append(status)
  footer = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
  save = Gtk.Button(label="Add button")
  close = Gtk.Button(label="Close")
  footer.append(save)
  footer.append(close)
  outer.append(footer)
  window.set_child(outer)
  editor = {"id": None, "source": None, "action": None, "icon": "", "menu_rows": [],
            "saved_rows": [], "icon_window": None, "active": False, "loading": False}
  def update_save(*_args):
    if editor["loading"]:
      return
    action_ready = (bool(call.get_text().strip()) if mode.get_active() == 1 else
                    bool(editor["action"]) if mode.get_active() == 2 else bool(editor["source"]))
    save.set_sensitive(editor["active"] and bool(label.get_text().strip()) and action_ready)
  def set_icon(token):
    editor["icon"] = token
    child = preview.get_first_child()
    while child is not None:
      preview.remove(child)
      child = preview.get_first_child()
    image = _toolbar_icon_image(token, window)
    preview.append(image if image is not None else Gtk.Label(label="Text only"))
    preview.set_tooltip_text(token or "No icon")
  def import_icon_file(path):
    if not editor["active"]:
      status.set_text("Click New button or choose a saved button first.")
      return False
    picker = icon_picker()
    return picker.toolbar_icon_widgets["import_path"](path)
  drop_target = _toolbar_add_icon_drop(icon_row, import_icon_file, status.set_text)
  if drop_target is not None:
    hint = Gtk.Label(label="Or drop an image onto the icon row.", xalign=0)
    hint.add_css_class("dim-label")
    right.insert_child_after(hint, icon_row)
  def menu_selected(_list, row):
    if not editor["active"] or editor["loading"]:
      return
    item = next((item for candidate, item in editor["menu_rows"] if candidate == row), None)
    if item is not None:
      if mode.get_active() == 2:
        old = editor["action"]
        editor["action"] = item["name"]
        if not label.get_text().strip() or label.get_text() == old:
          label.set_text(item["name"])
        source_label.set_text(item["name"] + " (" + _custom_key_binding_label(item["key"]) + ")")
        update_save()
        return
      old = editor["source"]
      editor["source"] = _toolbar_menu_source(item)
      if not label.get_text().strip() or (old and label.get_text() in (old["label"], old["label"].removesuffix("..."))):
        label.set_text(item["label"].removesuffix("..."))
      source_label.set_text(" > ".join(part for part in (item["path"], item["label"]) if part))
      update_save()
  def filter_menu(*_args):
    child = choices.get_first_child()
    while child is not None:
      choices.remove(child)
      child = choices.get_first_child()
    editor["menu_rows"] = []
    is_keybinding = mode.get_active() == 2
    items = _toolbar_keybinding_items(menu_search.get_text()) if is_keybinding else _toolbar_menu_items(menu_search.get_text())
    for item in items:
      row = Gtk.ListBoxRow()
      text = (item["name"] + " (" + _custom_key_binding_label(item["key"]) + ")" if is_keybinding else
              " > ".join(part for part in (item["path"], item["label"]) if part))
      item_label = Gtk.Label(label=text, xalign=0)
      item_label.set_wrap(True)
      item_label.set_max_width_chars(48)
      row.set_child(item_label)
      choices.append(row)
      editor["menu_rows"].append((row, item))
      if (item["name"] == editor["action"] if is_keybinding else _toolbar_menu_source(item) == editor["source"]):
        choices.select_row(row)
  def load(button=None, draft=False):
    editor["loading"] = True
    if editor["icon_window"] is not None:
      editor["icon_window"].destroy()
      editor["icon_window"] = None
    editor.update(id=None if button is None else button["id"], source=None if button is None else button.get("source"),
                  action=None if button is None else button.get("action"), active=button is not None or draft)
    heading.set_text("Editing: " + button["label"] if button is not None else "New button" if draft else "Choose a button")
    save.set_label("Add button" if button is None else "Update button")
    label.set_text("" if button is None else button["label"])
    tooltip.set_text("" if button is None else button["tooltip"])
    call.set_text("" if button is None else button.get("call", ""))
    show_label.set_active(button is None or button["show_label"])
    visible.set_active(button is None or button["visible"])
    mode.set_active({"call": 1, "keybinding": 2}.get(button["kind"], 0) if button else 0)
    set_icon("" if button is None else button["icon"])
    source = editor["source"]
    source_label.set_text(editor["action"] or ("" if source is None else source["path"] + " > " + source["label"]))
    filter_menu()
    for control in ("up", "down", "hide", "remove"):
      controls[control].set_sensitive(button is not None)
    # No GTK selection exists until an existing button is explicitly opened.
    saved.set_selection_mode(Gtk.SelectionMode.SINGLE if button is not None else Gtk.SelectionMode.NONE)
    row = next((row for row, ident in editor["saved_rows"] if ident == editor["id"]), None)
    if row is not None:
      saved.select_row(row)
    else:
      saved.unselect_all()
    right.set_sensitive(editor["active"])
    editor["loading"] = False
    update_save()
    status.set_text("New button: enter a label and choose an action." if draft else
        "Editing: " + button["label"] if button is not None else
        "Click New button to create one, or choose a saved button to edit.")
  def refresh_saved(selected_id=None):
    editor["loading"] = True
    child = saved.get_first_child()
    while child is not None:
      saved.remove(child)
      child = saved.get_first_child()
    editor["saved_rows"] = []
    selected = None
    for button in state["buttons"]:
      row = Gtk.ListBoxRow()
      item = Gtk.Label(label=button["label"] + (" (hidden)" if not button["visible"] else ""), xalign=0)
      item.set_wrap(True)
      item.set_max_width_chars(24)
      row.set_child(item)
      saved.append(row)
      editor["saved_rows"].append((row, button["id"]))
      if button["id"] == selected_id:
        selected = row
    saved.select_row(selected)
    editor["loading"] = False
  def saved_selected(_list, row):
    if editor["loading"]:
      return
    button_id = next((ident for candidate, ident in editor["saved_rows"] if candidate == row), None)
    if button_id is not None and button_id != editor["id"]:
      load(next(button for button in state["buttons"] if button["id"] == button_id))
  def commit(buttons, selected_id=None):
    try:
      _toolbar_commit(buttons)
      refresh_saved(selected_id)
      load(next((button for button in state["buttons"] if button["id"] == selected_id), None))
      status.set_text("Saved for this session and future starts.")
      return True
    except Exception as error:
      status.set_text("Could not save toolbar: " + str(error))
      return False
  def save_button(*_args):
    if not editor["active"]:
      return
    button = {"id": editor["id"] or os.urandom(12).hex(), "label": label.get_text(),
              "tooltip": tooltip.get_text(), "icon": editor["icon"], "show_label": show_label.get_active(),
              "visible": visible.get_active(), "kind": ("menu", "call", "keybinding")[mode.get_active()]}
    button.update({"call": call.get_text()} if button["kind"] == "call" else
                  {"action": editor["action"]} if button["kind"] == "keybinding" else {"source": editor["source"]})
    try:
      button = _toolbar_normalize_buttons([button])[0]
    except Exception as error:
      status.set_text(str(error))
      return
    buttons = [button if old["id"] == button["id"] else old for old in state["buttons"]]
    if editor["id"] is None:
      buttons.append(button)
    is_new = editor["id"] is None
    if commit(buttons, None if is_new else button["id"]):
      status.set_text(("Added '" if is_new else "Updated '") + button["label"] + "'. Saved for future starts.")
      problem = _toolbar_problem(button)
      if problem:
        status.set_text("Saved, but disabled in this session: " + problem)
  def change(operation):
    buttons = [dict(button) for button in state["buttons"]]
    index = next((i for i, button in enumerate(buttons) if button["id"] == editor["id"]), None)
    if index is None:
      return
    selected_id = editor["id"]
    if operation in ("up", "down"):
      other = index + (-1 if operation == "up" else 1)
      if not 0 <= other < len(buttons):
        return
      buttons[index], buttons[other] = buttons[other], buttons[index]
    elif operation == "hide":
      buttons[index]["visible"] = not buttons[index]["visible"]
    else:
      buttons.pop(index)
      selected_id = None
    if commit(buttons, selected_id) and selected_id is None:
      load()
  new = Gtk.Button(label="New button")
  new.connect("clicked", lambda *_: (load(draft=True), label.grab_focus()))
  left.append(new)
  order = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
  for key, text in (("up", "Up"), ("down", "Down"), ("hide", "Hide / Show"), ("remove", "Remove")):
    button = Gtk.Button(label=text)
    button.connect("clicked", lambda _widget, key=key: change(key))
    controls[key] = button
    order.append(button)
  left.append(order)
  def icon_picker(*_args):
    if not editor["active"]:
      return None
    previous = editor["icon_window"]
    if previous is not None and not previous.toolbar_icon_widgets["state"]["closed"]:
      previous.present()
      return previous
    editor["icon_window"] = _toolbar_choose_icon(window, set_icon, initial=editor["icon"],
        library_changed=lambda: set_icon(editor["icon"]))
    return editor["icon_window"]
  def closed(*_args):
    if editor["icon_window"] is not None:
      editor["icon_window"].destroy()
    state["window"] = None
    return False
  def mode_changed(*_args):
    is_call = mode.get_active() == 1
    pages.set_visible_child_name("call" if is_call else "menu")
    source_label.set_visible(not is_call)
    menu_search.set_placeholder_text("Find a keybinding action..." if mode.get_active() == 2 else "Find a menu action...")
    source_label.set_text("")
    filter_menu()
    update_save()
  mode.connect("changed", mode_changed)
  menu_search.connect("search-changed", filter_menu)
  choices.connect("row-selected", menu_selected)
  choices.connect("row-activated", menu_selected)
  # In blank/new mode rows are activatable but not selectable. While editing,
  # keyboard and mouse selection both load the highlighted button.
  saved.connect("row-activated", saved_selected)
  saved.connect("row-selected", saved_selected)
  label.connect("changed", update_save)
  call.connect("changed", update_save)
  choose_icon.connect("clicked", icon_picker)
  clear_icon.connect("clicked", lambda *_: set_icon(""))
  save.connect("clicked", save_button)
  close.connect("clicked", lambda *_: window.close())
  window.connect("close-request", closed)
  state["window"] = window
  window.toolbar_widgets = dict(controls, saved=saved, label=label, tooltip=tooltip, mode=mode, call=call,
      menu_search=menu_search, choices=choices, save=save, new=new, visible=visible, show_label=show_label,
      choose_icon=choose_icon, editor=editor, set_icon=set_icon, status=status, filter_menu=filter_menu,
      heading=heading, drop_target=drop_target, form=right)
  refresh_saved()
  load()
  if state["load_error"]:
    status.set_text(state["load_error"])
  window.set_focus(new)
  window.present()
  return window


def _install_custom_toolbar(menu):
  state = _toolbar_state()
  state["menu"] = menu
  toolbar = _coot_gui_module().coot_gui_api.main_toolbar()
  group = getattr(toolbar, "coot_trimmings_buttons", None)
  if group is None:
    group = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=4)
    toolbar.append(group)
    toolbar.coot_trimmings_buttons = group
  state["group"] = group
  _toolbar_render()


def build_custom_menu():
  """Create the top-level Custom menu with direct Python submenu builders."""
  coot_gui_module = _coot_gui_module()
  if coot_gui_module is None:
    return None
  if "Gio" not in globals() and hasattr(coot_gui_module, "Gio"):
    globals()["Gio"] = coot_gui_module.Gio

  menu = attach_module_menu_button("Custom")
  add_simple_coot_menu_menuitem(menu, "Custom Key Bindings...",
                               lambda _item: custom_key_bindings_dialog())
  add_simple_coot_menu_menuitem(menu, "Customize toolbar...",
                               lambda _item: custom_toolbar_dialog())
  submenu_display = Gio.Menu.new()
  submenu_colour = Gio.Menu.new()
  submenu_fit = Gio.Menu.new()
  submenu_renumber = Gio.Menu.new()
  submenu_settings = Gio.Menu.new()
  submenu_build = Gio.Menu.new()
  submenu_common_monomers = Gio.Menu.new()
  submenu_common_monomer_favorites = Gio.Menu.new()
  submenu_coordination_links = Gio.Menu.new()
  submenu_covalent_modifications = Gio.Menu.new()
  submenu_mutate = Gio.Menu.new()
  submenu_modify = Gio.Menu.new()
  submenu_maps = Gio.Menu.new()
  submenu_recent_models = Gio.Menu.new()
  submenu_recent_maps = Gio.Menu.new()

  menu.append_submenu("Display", submenu_display)
  if _supports_direct_user_defined_colouring():
    menu.append_submenu("User-defined colouring", submenu_colour)
  menu.append_submenu("Fit", submenu_fit)
  menu.append_submenu("Renumber", submenu_renumber)
  menu.append_submenu("Settings", submenu_settings)
  menu.append_submenu("Build", submenu_build)
  menu.append_submenu("Mutate", submenu_mutate)
  menu.append_submenu("Modify", submenu_modify)
  menu.append_submenu("Maps", submenu_maps)
  menu.append_submenu("Recent maps", submenu_recent_maps)
  menu.append_submenu("Recent models", submenu_recent_models)

  _build_custom_display_menu(submenu_display)
  if _supports_direct_user_defined_colouring():
    _build_custom_user_colour_menu(submenu_colour)
  _build_custom_fit_menu(submenu_fit)
  _build_custom_renumber_menu(submenu_renumber)
  _build_custom_settings_menu(submenu_settings)
  _build_custom_build_menu(
    submenu_build,
    submenu_common_monomers,
    submenu_common_monomer_favorites,
    submenu_coordination_links,
    submenu_covalent_modifications,
  )
  _build_custom_mutate_menu(submenu_mutate)
  _build_custom_modify_menu(submenu_modify)
  _build_custom_maps_menu(submenu_maps)
  _install_recent_files_tracking(submenu_recent_models, submenu_recent_maps)
  try:
    _install_custom_menu_search(menu)
  except Exception as error:
    print("Custom menu search unavailable; keeping the standard menu:", error)
  try:
    _install_custom_toolbar(menu)
  except Exception as error:
    print("Custom toolbar unavailable; keeping Coot's toolbar unchanged:", error)
  return menu

if _coot_gui_module() is not None:
  try:
    build_custom_menu()
  except Exception:
    print("coot_trimmings direct menu build failed")
    traceback.print_exc()
  try:
    _install_trackpad_native_close_guard()
  except Exception as error:
    print("Native Rotate/Translate close guard unavailable:", error)

_initialize_trackpad_preferences()

# Keep this as the final startup action in the file.  It covers settings Coot
# may have changed while the rest of this startup script was being evaluated;
# the GLib callbacks above then cover state restoration after script completion.
_apply_persistent_display_defaults()
