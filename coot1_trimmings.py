"""Oli's Coot 1.x customizations.

Copy this file into Coot's startup script directory and restart Coot.
The main user-editable knobs live near the top of the file under
"User-editable settings".
"""

import math
import base64
import hashlib
import json
import os
import re
import shutil
import struct
import sys
import tempfile
import traceback
import warnings
from datetime import datetime, timezone
from functools import lru_cache

STARTUP_DIAGNOSTIC_MESSAGES = []
_STARTUP_DIAGNOSTIC_DIALOGS_SHOWN = set()
_GTK_REPOSITORY_WARNINGS_SHOWN = set()


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
  coot_gui_module = _coot_gui_module()
  if coot_gui_module is not None and hasattr(coot_gui_module, module_name):
    return getattr(coot_gui_module, module_name)
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
STARTUP_VARIABLE_BOND_THICKNESS = 7
STARTUP_DEFAULT_NEW_ATOM_B = 50.0
STARTUP_ADD_TERMINAL_POST_REFINE = 1
STARTUP_TERMINAL_RIGID_BODY_REFINE = 0
STARTUP_MUTATE_AUTO_FIT_POST_REFINE = 1
STARTUP_SYMMETRY_SIZE = 30
STARTUP_NOMENCLATURE_ERRORS_MODE = "ignore"
STARTUP_REORIENTING_NEXT_RESIDUE_MODE = 1
STARTUP_ROTATION_CENTRE_CROSSHAIRS_SCALE = 0.15
STARTUP_ROTATION_CENTRE_CROSSHAIRS_COLOUR = (1.0, 1.0, 1.0, 1.0)

# User-facing helper defaults.
DEFAULT_NEW_HELIX_CHAIN_ID = "A"
HIGH_CONTRAST_BOND_THICKNESS = 2
COMMON_MONOMER_FAVORITES_FILENAME = "coot_trimmings_favorites.json"
COMMON_MONOMER_FAVORITE_CIF_PREFIX = "coot_trimmings_favorite_"
WATER_REFINE_MAX_DISPLACEMENT = 1.0
RESIDUE_ANNOTATION_DEFAULT_AUTHOR = (
  os.environ.get("USER")
  or os.environ.get("USERNAME")
  or ""
)

# Map restyling defaults for the EM helper.
EM_REFINED_MAP_COLOUR = (0.10, 0.57, 0.95)
EM_REFINED_MAP_CONTOUR_SIGMA = 2.3
EM_TARGET_PIXEL_SIZE = 0.5
EM_RESAMPLE_CONFIRM_MAX_GRID_DIMENSION = 1024
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
MAP_SURFACE_DISPLAY_STATE = {}
MAP_SURFACE_OPACITY_STATE = {}
MAP_GLOBAL_VIEW_SETTINGS = {}
MAP_LOCAL_APPEARANCE_STATE = {}
SMART_COPY_TEMPLATE_IMOL = None
SMART_COPY_SOURCE_CENTRE = None
SMART_COPY_RESIDUE_NAME = None
MODEL_AMBIENT_LIGHTING_ENABLED = False
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
      imol_backwards = coot.copy_molecule(imol)
      if valid_model_molecule_qm(imol_backwards):
        temp_imols.add(imol_backwards)
      loop_len = abs(start_resno - stop_resno) + 1
      imol_both = coot.copy_molecule(imol) if loop_len >= 6 else None
      if imol_both is not None and valid_model_molecule_qm(imol_both):
        temp_imols.add(imol_both)

      atom_selection = "//" + chain_id + "/" + str(min(start_resno, stop_resno) - 1) + \
                       "-" + str(max(start_resno, stop_resno) + 1)
      imol_fragment_backup = coot.new_molecule_by_atom_selection(imol, atom_selection)
      if valid_model_molecule_qm(imol_fragment_backup):
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
        coot.set_refinement_immediate_replacement(1)
        coot.refine_zone(imol_both, chain_id, stop_resno1 - loop_len//3, start_resno2 + loop_len//3, "")
        coot.accept_regularizement()
        coot.set_refinement_immediate_replacement(immediate_refinement_mode)
        result_c = gap.low_density_average(imol_map, imol_both, chain_id, start_resno, stop_resno)
        loop_list.append([imol_both, result_c])

      i = 0
      while i < (len(loop_list) - 1):
        j = i + 1
        while j < len(loop_list):
          max_score = max(loop_list[i][1], loop_list[j][1])
          if max_score != 0 and (min(loop_list[i][1], loop_list[j][1]) / max_score > 0.90):
            temp_imols.discard(loop_list[j][0])
            coot.close_molecule(loop_list[j][0])
            loop_list.pop(j)
          j += 1
        i += 1

      if result_a > result_b:
        temp_imols.discard(imol_backwards)
        coot.close_molecule(imol_backwards)
      else:
        coot.replace_fragment(imol, imol_backwards, atom_selection)
        temp_imols.discard(imol_backwards)
        coot.close_molecule(imol_backwards)
    else:
      gap.fit_gap_generic(imol, chain_id, start_resno, stop_resno, sequence)
  finally:
    for temp_imol in list(temp_imols):
      if valid_model_molecule_qm(temp_imol):
        coot.close_molecule(temp_imol)
    coot.set_refine_ramachandran_angles(rama_status)
    if backup_mode == 1:
      coot.turn_on_backup(imol)

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


# ============================================================================
# Legacy / disabled custom-colouring support
# ============================================================================
#
# The Colour menu is currently disabled in Coot 1.1 because the user-defined
# colouring APIs remain awkward and fragile. The code below is kept for
# reference and for a few programmatic helpers, but it is intentionally kept
# away from the main startup settings and general utility code.
#

USER_DEFINED_COLOUR_TABLE_SIZE = 60
USER_DEFINED_COLOUR_TABLE_BASE = 60
USER_DEFINED_COLOUR_TABLE_READY = False


def _make_default_user_defined_colours():
  base_colours = [[0.72, 0.72, 0.72] for _ in range(USER_DEFINED_COLOUR_TABLE_SIZE)]
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
  }
  for colour_index, colour in explicit_colours.items():
    if 0 <= colour_index < USER_DEFINED_COLOUR_TABLE_SIZE:
      base_colours[colour_index] = colour

  colours = []
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


def ensure_user_defined_colour_table():
  global USER_DEFINED_COLOUR_TABLE_READY
  if USER_DEFINED_COLOUR_TABLE_READY:
    return
  if "make_alphafold_colours" in globals():
    colours = make_alphafold_colours()
  else:
    colours = _make_default_user_defined_colours()
  coot.set_user_defined_colours_py(colours)
  USER_DEFINED_COLOUR_TABLE_READY = True


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


_coot_graphics_to_bonds_representation = graphics_to_bonds_representation
_coot_graphics_to_rainbow_representation = graphics_to_rainbow_representation
_coot_graphics_to_b_factor_representation = graphics_to_b_factor_representation
_coot_graphics_to_ca_plus_ligands_representation = graphics_to_ca_plus_ligands_representation
_coot_graphics_to_ca_plus_ligands_and_sidechains_representation = (
  graphics_to_ca_plus_ligands_and_sidechains_representation
)
_coot_graphics_to_ca_plus_ligands_sec_struct_representation = (
  graphics_to_ca_plus_ligands_sec_struct_representation
)


def _clear_user_defined_colours_for_standard_representation(mol_id):
  if mol_id is None or mol_id < 0:
    return
  clear_user_defined_atom_colours(mol_id)


def graphics_to_bonds_representation(mol_id):
  _clear_user_defined_colours_for_standard_representation(mol_id)
  return _coot_graphics_to_bonds_representation(mol_id)


def graphics_to_rainbow_representation(mol_id):
  _clear_user_defined_colours_for_standard_representation(mol_id)
  return _coot_graphics_to_rainbow_representation(mol_id)


def graphics_to_b_factor_representation(mol_id):
  _clear_user_defined_colours_for_standard_representation(mol_id)
  return _coot_graphics_to_b_factor_representation(mol_id)


def graphics_to_ca_plus_ligands_representation(mol_id):
  _clear_user_defined_colours_for_standard_representation(mol_id)
  return _coot_graphics_to_ca_plus_ligands_representation(mol_id)


def graphics_to_ca_plus_ligands_and_sidechains_representation(mol_id):
  _clear_user_defined_colours_for_standard_representation(mol_id)
  return _coot_graphics_to_ca_plus_ligands_and_sidechains_representation(mol_id)


def graphics_to_ca_plus_ligands_sec_struct_representation(mol_id):
  _clear_user_defined_colours_for_standard_representation(mol_id)
  return _coot_graphics_to_ca_plus_ligands_sec_struct_representation(mol_id)


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
  "A", "C", "G", "U", "T", "DA", "DC", "DG", "DT",
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
  "Weak Sidechains",
  "Weak Backbone",
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
    # Older Coot builds may not expose the variable-thickness API at all.
    set_use_variable_bond_thickness(STARTUP_USE_VARIABLE_BOND_THICKNESS)
    set_default_bond_thickness(STARTUP_VARIABLE_BOND_THICKNESS)
  except NameError:
    print("Your coot is a bit old... consider upgrading...")
  set_default_temperature_factor_for_new_atoms(STARTUP_DEFAULT_NEW_ATOM_B)
  set_add_terminal_residue_do_post_refine(STARTUP_ADD_TERMINAL_POST_REFINE)
  set_terminal_residue_do_rigid_body_refine(STARTUP_TERMINAL_RIGID_BODY_REFINE)
  set_mutate_auto_fit_do_post_refine(STARTUP_MUTATE_AUTO_FIT_POST_REFINE)
  set_symmetry_size(STARTUP_SYMMETRY_SIZE)
  set_nomenclature_errors_on_read(STARTUP_NOMENCLATURE_ERRORS_MODE)
  set_reorienting_next_residue_mode(STARTUP_REORIENTING_NEXT_RESIDUE_MODE)
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


def _segment_range_for_residue(mol_id, chain_id, resno):
  """Return the contiguous polymer segment containing a residue number."""
  for _mol_id, seg_chain_id, seg_start, seg_end in segment_list_chain(mol_id, chain_id):
    if seg_chain_id == chain_id and seg_start <= resno <= seg_end:
      return (seg_start, seg_end)
  return None


def _active_segment_context_or_status():
  """Return active atom/segment context for segment-level helpers."""
  atom_context = _active_atom_context_or_status()
  if not atom_context:
    return None
  segment_range = _segment_range_for_residue(
    atom_context["mol_id"],
    atom_context["chain_id"],
    atom_context["resno"],
  )
  if segment_range is None:
    add_status_bar_text("No active segment")
    return None
  atom_context["segment_start"] = segment_range[0]
  atom_context["segment_end"] = segment_range[1]
  return atom_context


def _refresh_extra_restraints_display(mol_id):
  """Rebuild the displayed extra-restraints state after editing/renumbering."""
  delete_all_extra_restraints(mol_id)
  set_show_extra_restraints(mol_id, 0)
  set_show_extra_restraints(mol_id, 1)


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
  chain_ids_to_search = [chain_id] if chain_id and chain_id in chain_ids(mol_id) else chain_ids(mol_id)
  for current_chain_id in chain_ids_to_search:
    for serial_number in range(chain_n_residues(current_chain_id, mol_id)):
      resno = seqnum_from_serial_number(mol_id, current_chain_id, serial_number)
      ins_code = insertion_code_from_serial_number(mol_id, current_chain_id, serial_number)
      if not _residue_is_polymer(mol_id, current_chain_id, resno, ins_code):
        continue
      residue_centre = residue_centre_py(mol_id, current_chain_id, resno, ins_code)
      if not isinstance(residue_centre, (list, tuple)) or len(residue_centre) != 3:
        continue
      distance_sq = _distance_sq(centre, residue_centre)
      if distance_sq <= best_distance_sq:
        best_distance_sq = distance_sq
        best_residue = (current_chain_id, resno, ins_code)
  return best_residue


def _navigation_reference_residue():
  """Resolve the residue that navigation should use as its starting point."""
  residue = active_residue()
  if residue and _residue_is_polymer(residue[0], residue[1], residue[2], residue[3]):
    return {
      "mol_id": residue[0],
      "chain_id": residue[1],
      "resno": residue[2],
      "ins_code": residue[3],
      "source": "active",
    }

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
      return {
        "mol_id": mol_id,
        "chain_id": chain_id,
        "resno": resno,
        "ins_code": ins_code,
        "source": "stored",
        "serial_number": serial_number,
      }
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
  }
  return 1


def _rotation_centre_xyz():
  return [rotation_centre_position(axis) for axis in (0, 1, 2)]


def _distance_sq(point_1, point_2):
  dx = point_1[0] - point_2[0]
  dy = point_1[1] - point_2[1]
  dz = point_1[2] - point_2[2]
  return dx*dx + dy*dy + dz*dz


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
  for imol in model_molecule_list():
    if residue_exists_qm(imol, chain_id, resno, ins_code):
      return imol
  return -1


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
  """Return a best-effort model molecule id from a Coot click spec.

  Coot 1.x click payloads often contain correct chain/residue information but
  unreliable molecule ids. Prefer a valid model id when present, otherwise fall
  back to a residue lookup and finally the active residue.
  """
  valid_model_mols = model_molecule_list()
  if not isinstance(click_spec, (list, tuple)):
    residue = active_residue()
    return residue[0] if residue and residue[0] in valid_model_mols else -1
  first_imol = _coot_int_arg(_click_spec_field(click_spec, 0, 0, -1))
  second_imol = _coot_int_arg(_click_spec_field(click_spec, 1, 0, -1))
  if isinstance(first_imol, int) and first_imol in valid_model_mols:
    return first_imol
  if len(click_spec) >= 7 and isinstance(second_imol, int) and second_imol in valid_model_mols:
    return second_imol
  chain_id = _click_spec_chain_id(click_spec)
  resno = _click_spec_res_no(click_spec)
  ins_code = _click_spec_ins_code(click_spec)
  imol = _find_model_molecule_for_click_spec(chain_id, resno, ins_code)
  if imol in valid_model_mols:
    return imol
  residue = active_residue()
  if residue and residue[0] in valid_model_mols:
    return residue[0]
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

  def _populate_entry_from_file_chooser(parent_window, target_entry, chooser_title):
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
      initial_directory = COMMON_MONOMER_LAST_BROWSED_DIRECTORY
      if not initial_directory or not os.path.isdir(initial_directory):
        return None
      if gio_module is None:
        gio_module = _gtk_gio_module()
      if gio_module is None or not hasattr(gio_module, "File"):
        return None
      initial_folder = gio_module.File.new_for_path(initial_directory)
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
        chooser.set_default_filter(cif_filter)

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

    info_dialog("Browse is unavailable in this Gtk build; please paste the CIF path manually.")
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
  ):
    """GTK4-safe one-entry prompt with a Browse button for save/export paths."""
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
      lambda *_args: _populate_entry_from_save_file_chooser(window, entry, entry_1_browse_title),
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

  def generic_author_note_entry(
    function_label,
    residue_label,
    author_default_text,
    note_default_text,
    go_button_label,
    handle_go_function,
  ):
    """Small GTK4 prompt for adding an authored multiline residue note."""
    window = Gtk.Window()
    window.set_title("Coot")
    window.set_default_size(520, 360)

    vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=8)
    hbox_buttons = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
    title_label = Gtk.Label(label=function_label)
    residue_info_label = Gtk.Label(label=residue_label)
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
    hbox_buttons.set_margin_start(12)
    hbox_buttons.set_margin_end(12)
    hbox_buttons.set_margin_bottom(12)

    if isinstance(author_default_text, str):
      author_entry.set_text(author_default_text)
    if isinstance(note_default_text, str):
      note_buffer.set_text(note_default_text)

    note_view.set_wrap_mode(Gtk.WrapMode.WORD_CHAR)
    note_scrolled.set_child(note_view)

    def close_window(*_args):
      window.destroy()
      return False

    def submit(*_args):
      start_iter = note_buffer.get_start_iter()
      end_iter = note_buffer.get_end_iter()
      status = handle_go_function(
        author_entry.get_text(),
        note_buffer.get_text(start_iter, end_iter, False),
      )
      if status not in (0, False):
        window.destroy()
      return False

    cancel_button.connect("clicked", close_window)
    go_button.connect("clicked", submit)
    author_entry.connect("activate", lambda *_args: note_view.grab_focus())

    vbox.append(title_label)
    vbox.append(residue_info_label)
    vbox.append(author_label)
    vbox.append(author_entry)
    vbox.append(note_label)
    vbox.append(note_scrolled)
    hbox_buttons.append(cancel_button)
    hbox_buttons.append(go_button)
    vbox.append(hbox_buttons)
    window.set_child(vbox)
    window.present()
    if author_entry.get_text().strip():
      note_view.grab_focus()
    else:
      author_entry.grab_focus()

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

  def categorized_interesting_things_gui(dialog_name, categorized_thing_lists, refresh_function=None):
    """Show grouped clickable jump targets in a single scrollable dialog."""
    window = Gtk.Window()
    window.set_title("Coot")
    window.set_default_size(560, 520)

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
      window.destroy()

    def jump_to_entry(entry, category_name=None):
      _activate_interesting_entry(entry, category_name)

    # Keep a flat ordered list for Prev/Next while still rendering grouped rows.
    # `None` means the user has not navigated the list yet.
    flat_entry_specs = []
    current_index = {"value": None}
    selected_button = {"value": None}

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

    def render_categories(current_dialog_name, current_categorized_thing_lists):
      label.set_text(current_dialog_name)
      current_index["value"] = None
      set_selected_button(None)
      flat_entry_specs.clear()

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

        header = Gtk.Label(label=f"{category_name} ({len(thing_list)})")
        header.set_xalign(0.0)
        header.set_margin_top(6)
        header.set_margin_bottom(2)
        inside_vbox.append(header)

        for entry in thing_list:
          if len(entry) < 4:
            continue
          button = Gtk.Button(label=str(entry[0]))
          entry_index = len(flat_entry_specs)
          flat_entry_specs.append((category_name, entry, button))
          button.connect(
            "clicked",
            lambda _button, idx=entry_index: activate_entry(idx),
          )
          inside_vbox.append(button)

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
      refreshed_payload = refresh_function()
      if not refreshed_payload:
        return None
      refreshed_dialog_name, refreshed_categorized_thing_lists = refreshed_payload
      render_categories(refreshed_dialog_name, refreshed_categorized_thing_lists)

    prev_button.connect("clicked", lambda *_args: activate_prev())
    next_button.connect("clicked", lambda *_args: activate_next())
    refresh_button.connect("clicked", refresh_categories)
    close_button.connect("clicked", close_window)
    refresh_button.set_sensitive(callable(refresh_function))

    render_categories(dialog_name, categorized_thing_lists)

    vbox.append(label)
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
  ):
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
  ):
    info_dialog(_gui_unavailable_message(function_label))

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

  def categorized_interesting_things_gui(dialog_name, categorized_thing_lists, refresh_function=None):
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

  def wrapped_thunk():
    try:
      return thunk()
    except Exception:
      print("coot_trimmings keybinding failure:", name, "key=", key_value)
      traceback.print_exc()
      return None

  REGISTERED_KEYBINDING_CALLBACKS.append(wrapped_thunk)
  coot.add_key_binding_gtk4_py(key_value, ctrl_key, wrapped_thunk, name)

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
lambda: clear_pending_picks())

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
      default_bond_thickness = (
        STARTUP_VARIABLE_BOND_THICKNESS
        if STARTUP_USE_VARIABLE_BOND_THICKNESS
        else STARTUP_DEFAULT_BOND_THICKNESS
      )
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
    copied_imol = new_molecule_by_atom_selection(mol_id, atom_selection)
  finally:
    set_recentre_on_read_pdb(recentre_state)
    _restore_view_state(view_state)
  if copied_imol == -1:
    add_status_bar_text("Unable to copy active residue")
    return None
  _close_smart_copy_template()
  _set_smart_copy_template(copied_imol)
  SMART_COPY_SOURCE_CENTRE = residue_centre_py(mol_id, ch_id, resno, ins_code)
  SMART_COPY_RESIDUE_NAME = residue_name(mol_id, ch_id, resno, ins_code)
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
  molecule_label = molecule_name_stub_py(target_mol_id, 0) or molecule_name(target_mol_id) or ""
  molecule_label = os.path.basename(str(molecule_label).strip()) or f"molecule_{target_mol_id}"
  return f"Molecule #{target_mol_id}, {molecule_label}"


def _merge_molecules_explicitly_failed(merge_result):
  if not isinstance(merge_result, (list, tuple)) or not merge_result:
    return False
  return merge_result[0] in (0, False, None)


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
    "Pasted {resname} as a new molecule ({target})".format(
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
    merge_result = merge_molecules_py([paste_imol], target_mol_id)
    if _merge_molecules_explicitly_failed(merge_result):
      add_status_bar_text("Unable to merge pasted residue into active model")
      return None
  finally:
    if paste_imol in molecule_number_list():
      close_molecule(paste_imol)
  add_status_bar_text(
    "Pasted {resname} at pointer in {target}".format(
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


def _residue_atom_records_and_xyz(imol, chain_id, resno, ins_code):
  """Parse a residue once into both atom records and an atom-name lookup map."""
  atom_info = residue_info_py(imol, chain_id, resno, ins_code)
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


def _collect_odd_residue_dialog_data():
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

  for chain_index, chain_id in enumerate(chain_ids(mol_id)):
    for serial_number in range(chain_n_residues(chain_id, mol_id)):
      residue_name_here = resname_from_serial_number(mol_id, chain_id, serial_number)
      resno = seqnum_from_serial_number(mol_id, chain_id, serial_number)
      ins_code = insertion_code_from_serial_number(mol_id, chain_id, serial_number)

      atom_records, atom_xyz = _residue_atom_records_and_xyz(mol_id, chain_id, resno, ins_code)
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
          navigation_metadata if use_navigation else None,
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

      if residue_name_here == "GLY":
        continue

      support_fraction, supported_points, total_points = _backbone_density_support(map_id, atom_xyz, peak_threshold)
      if total_points == 0 or support_fraction < EMRINGER_HELPER_MIN_BACKBONE_SUPPORT_FRACTION:
        add_result(
          "Weak Backbone",
          (support_fraction, chain_index, serial_number),
          f"{residue_id} {residue_name_here}: backbone support {supported_points}/{total_points} "
          f"({support_fraction:.0%}) below threshold",
          use_navigation=True,
        )
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

  sorted_gui_categories, sorted_detail_categories = _sorted_categorized_outputs(
    categorized_entries,
    ODD_RESIDUE_CATEGORY_ORDER,
  )
  numbered_gui_categories = _numbered_gui_categories(sorted_gui_categories)
  total_hits = sum(len(entries) for entries in categorized_entries.values())
  return {
    "mol_id": mol_id,
    "map_id": map_id,
    "peak_threshold": peak_threshold,
    "gui_categories": numbered_gui_categories,
    "detail_categories": sorted_detail_categories,
    "title": _odd_residue_dialog_title(mol_id, map_id),
    "total_hits": total_hits,
  }


def find_odd_residues():
  """Find odd residues in the active molecule using density-based heuristics."""
  dialog_data = _collect_odd_residue_dialog_data()
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
    return refreshed_data["title"], refreshed_data["gui_categories"]

  categorized_interesting_things_gui(
    dialog_data["title"],
    dialog_data["gui_categories"],
    refresh_function=refresh_payload,
  )
  _log_odd_residue_results(
    dialog_data["mol_id"],
    dialog_data["map_id"],
    dialog_data["peak_threshold"],
    dialog_data["detail_categories"],
  )
  return dialog_data["total_hits"]


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


def _annotation_safe_int(value, default=None):
  try:
    return int(str(value).strip())
  except Exception:
    return default


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


def _annotation_prepare_note_fields(author_text, note_text):
  normalized_author = _annotation_text_or_empty(author_text).strip() or _annotation_default_author()
  normalized_note = _annotation_text_or_empty(note_text).strip()
  if not normalized_note:
    info_dialog("Please enter a note before saving the residue annotation.")
    return None
  return normalized_author, normalized_note


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


def _annotation_normalize_entry(raw_entry, fallback_id=None):
  auth_asym_id = _annotation_text_or_empty(raw_entry.get("auth_asym_id")).strip()
  auth_seq_id = _annotation_safe_int(raw_entry.get("auth_seq_id"), default=None)
  if not auth_asym_id or auth_seq_id is None:
    return None
  entry_id = _annotation_safe_int(raw_entry.get("id"), default=fallback_id)
  if entry_id is None:
    return None
  normalized_entry = {
    "id": entry_id,
    "auth_asym_id": auth_asym_id,
    "auth_seq_id": auth_seq_id,
    "pdbx_PDB_ins_code": _annotation_text_or_empty(raw_entry.get("pdbx_PDB_ins_code")).strip(),
    "label_comp_id": _annotation_text_or_empty(raw_entry.get("label_comp_id")).strip() or "?",
    "author": _annotation_text_or_empty(raw_entry.get("author")).strip() or _annotation_default_author(),
    "modified_utc": _annotation_text_or_empty(raw_entry.get("modified_utc")).strip() or _annotation_now_utc(),
    "note": _annotation_text_or_empty(raw_entry.get("note")).strip(),
  }
  if not normalized_entry["note"]:
    return None
  return normalized_entry


def _annotation_group_sort_key(group):
  first_entry = group["entries"][0]
  return (group["key"][0], group["key"][1], group["key"][2], first_entry.get("id", 0))


def _annotation_entry_sort_key(entry):
  return (entry.get("modified_utc", ""), entry.get("id", 0))


def _annotation_groups_for_molecule(mol_id):
  grouped = {}
  for entry in _annotation_entries_for_molecule(mol_id):
    group_key = _annotation_residue_key(
      entry["auth_asym_id"],
      entry["auth_seq_id"],
      entry.get("pdbx_PDB_ins_code", ""),
      entry.get("label_comp_id", ""),
    )
    grouped.setdefault(group_key, []).append(entry)

  groups = []
  for group_key, entries in grouped.items():
    entries.sort(key=_annotation_entry_sort_key)
    chain_id, resno, ins_code, residue_name_here = group_key
    groups.append({
      "key": group_key,
      "dialog_label": _odd_residue_dialog_label(chain_id, resno, ins_code, residue_name_here),
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
    lines.append(f"Author: {_annotation_text_or_empty(entry.get('author')) or 'Unknown'}")
    lines.append(f"Timestamp: {_annotation_text_or_empty(entry.get('modified_utc')) or 'Unknown'}")
    lines.append("")
    lines.append(_annotation_text_or_empty(entry.get("note")) or "")
    if entry_index != len(group["entries"]):
      lines.extend(["", "-" * 40, ""])
  return "\n".join(lines)


def _annotation_group_preview_text(group, max_length=56):
  latest_note = _annotation_text_or_empty(group["entries"][-1].get("note")).splitlines()
  preview_text = latest_note[0].strip() if latest_note else ""
  if len(preview_text) > max_length:
    preview_text = preview_text[:max_length - 1].rstrip() + "..."
  note_count = len(group["entries"])
  suffix = "note" if note_count == 1 else "notes"
  if preview_text:
    return f"{group['dialog_label']} ({note_count} {suffix}) - {preview_text}"
  return f"{group['dialog_label']} ({note_count} {suffix})"


def _annotation_mmcif_rows(mol_id):
  rows = []
  for entry in sorted(_annotation_entries_for_molecule(mol_id), key=_annotation_entry_sort_key):
    rows.append(
      [
        _annotation_text_or_empty(entry.get("id")) or "?",
        _annotation_text_or_empty(entry.get("auth_asym_id")) or "?",
        str(entry.get("auth_seq_id", "")),
        _annotation_text_or_empty(entry.get("pdbx_PDB_ins_code")) or "?",
        _annotation_text_or_empty(entry.get("label_comp_id")) or "?",
        _annotation_encode_text(entry.get("author")),
        _annotation_text_or_empty(entry.get("modified_utc")) or "?",
        _annotation_encode_text(entry.get("note")),
      ]
    )
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
    "| Residue | Author | Timestamp | Note |",
    "| --- | --- | --- | --- |",
  ]
  for group in _annotation_groups_for_molecule(mol_id):
    residue_label = _annotation_markdown_cell(group["dialog_label"])
    for entry in group["entries"]:
      lines.append(
        "| {0} | {1} | {2} | {3} |".format(
          residue_label,
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
    RESIDUE_ANNOTATION_CATEGORY_PREFIX + "auth_asym_id",
    RESIDUE_ANNOTATION_CATEGORY_PREFIX + "auth_seq_id",
    RESIDUE_ANNOTATION_CATEGORY_PREFIX + "pdbx_PDB_ins_code",
    RESIDUE_ANNOTATION_CATEGORY_PREFIX + "label_comp_id",
    RESIDUE_ANNOTATION_CATEGORY_PREFIX + "modified_utc",
  ]
  for required_tag in required_tags:
    if required_tag not in tag_indices:
      raise RuntimeError(f"Embedded annotation loop is missing required tag {required_tag}.")

  author_tag = None
  for candidate_tag in (
    RESIDUE_ANNOTATION_CATEGORY_PREFIX + "author_b64",
    RESIDUE_ANNOTATION_CATEGORY_PREFIX + "author",
  ):
    if candidate_tag in tag_indices:
      author_tag = candidate_tag
      break
  note_tag = None
  for candidate_tag in (
    RESIDUE_ANNOTATION_CATEGORY_PREFIX + "note_b64",
    RESIDUE_ANNOTATION_CATEGORY_PREFIX + "note",
  ):
    if candidate_tag in tag_indices:
      note_tag = candidate_tag
      break
  if author_tag is None or note_tag is None:
    raise RuntimeError("Embedded annotation loop is missing author or note text fields.")

  entries = []
  for row_index in range(loop.length()):
    normalized_entry = _annotation_normalize_entry(
      {
        "id": _annotation_text_or_empty(loop[row_index, tag_indices[RESIDUE_ANNOTATION_CATEGORY_PREFIX + "id"]]),
        "auth_asym_id": _annotation_text_or_empty(loop[row_index, tag_indices[RESIDUE_ANNOTATION_CATEGORY_PREFIX + "auth_asym_id"]]),
        "auth_seq_id": _annotation_text_or_empty(loop[row_index, tag_indices[RESIDUE_ANNOTATION_CATEGORY_PREFIX + "auth_seq_id"]]),
        "pdbx_PDB_ins_code": _annotation_text_or_empty(loop[row_index, tag_indices[RESIDUE_ANNOTATION_CATEGORY_PREFIX + "pdbx_PDB_ins_code"]]),
        "label_comp_id": _annotation_text_or_empty(loop[row_index, tag_indices[RESIDUE_ANNOTATION_CATEGORY_PREFIX + "label_comp_id"]]),
        "author": _annotation_decode_text(loop[row_index, tag_indices[author_tag]]),
        "modified_utc": _annotation_text_or_empty(loop[row_index, tag_indices[RESIDUE_ANNOTATION_CATEGORY_PREFIX + "modified_utc"]]),
        "note": _annotation_decode_text(loop[row_index, tag_indices[note_tag]]),
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
  block.set_pair(RESIDUE_ANNOTATION_META_VERSION_TAG, "1")
  block.set_pair(RESIDUE_ANNOTATION_META_PROGRAM_TAG, "coot-trimmings")
  if entries:
    loop = block.init_loop(
      RESIDUE_ANNOTATION_CATEGORY_PREFIX,
      [
        "id",
        "auth_asym_id",
        "auth_seq_id",
        "pdbx_PDB_ins_code",
        "label_comp_id",
        "author_b64",
        "modified_utc",
        "note_b64",
      ],
    )
    for row in entries:
      loop.add_row(row)


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
    entity_table.append_row([record["entity_id"], record["entity_type"]])
    struct_asym_table.append_row([record["asym_id"], record["entity_id"]])
    if record["kind"] == "polymer" and entity_poly_table is not None and entity_poly_seq_table is not None:
      entity_poly_table.append_row([
        record["entity_id"],
        record["polymer_type"] or "other",
        record["asym_id"],
        record["one_letter_sequence"] or "?",
      ])
      for index, mon_id in enumerate(record["sequence"], start=1):
        entity_poly_seq_table.append_row([record["entity_id"], str(index), mon_id or "?"])
    elif record["kind"] == "nonpolymer" and nonpoly_table is not None:
      nonpoly_table.append_row([record["entity_id"], record["comp_id"] or "?"])

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

  for row_index in range(len(group_column)):
    auth_chain_id = auth_chain_column[row_index]
    ins_code = ins_code_column[row_index]
    residue_key = (
      auth_chain_id,
      auth_seq_column[row_index],
      "" if ins_code in (".", "?") else ins_code,
      comp_id_column[row_index],
    )
    nonpoly_record = nonpoly_residue_map.get(residue_key)
    if nonpoly_record is not None:
      chain_column[row_index] = nonpoly_record["asym_id"]
      entity_column[row_index] = nonpoly_record["entity_id"]
      label_seq_column[row_index] = "."
      continue

    polymer_record = polymer_residue_map.get(auth_chain_id)
    if polymer_record is None:
      label_seq_column[row_index] = "."
      continue

    chain_column[row_index] = polymer_record["asym_id"]
    entity_column[row_index] = polymer_record["entity_id"]
    label_seq_column[row_index] = polymer_record["residue_key_to_label_seq"].get(residue_key, ".")


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
  return document


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
    export_document = _annotation_export_document_from_coot_cif(source_path)
    _annotation_write_entries_to_mmcif_block(export_document.sole_block(), _annotation_mmcif_rows(mol_id))
    export_handle, export_path = tempfile.mkstemp(suffix=".cif", dir=output_directory)
    os.close(export_handle)
    export_handle = None
    export_document.write_file(export_path)
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


def _append_annotation_for_target(target, author_text, note_text, refresh_function=None):
  global RESIDUE_ANNOTATION_LAST_AUTHOR
  if not target:
    return 0
  note_fields = _annotation_prepare_note_fields(author_text, note_text)
  if note_fields is None:
    return 0
  normalized_author, normalized_note = note_fields
  RESIDUE_ANNOTATION_LAST_AUTHOR = normalized_author

  mol_id = target["mol_id"]
  annotation_entry = _new_annotation_entry(
    target,
    normalized_author,
    normalized_note,
    _annotation_next_id(mol_id),
  )
  _annotation_entries_for_molecule(mol_id).append(annotation_entry)
  residue_label = _odd_residue_dialog_label(
    target["chain_id"],
    target["resno"],
    target["ins_code"],
    target["residue_name"],
  )
  add_status_bar_text(f"Added residue annotation for {residue_label}")
  if callable(refresh_function):
    refresh_function(
      _annotation_residue_key(
        target["chain_id"],
        target["resno"],
        target["ins_code"],
        target["residue_name"],
      )
    )
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
  generic_author_note_entry(
    "Add residue annotation",
    residue_label,
    _annotation_default_author(),
    "",
    "Add note",
    lambda author_text, note_text: _append_annotation_for_target(
      target,
      author_text,
      note_text,
      refresh_function,
    ),
  )
  return 1


def _annotation_note_choice_label(entry, entry_index, max_length=48):
  author_text = _annotation_text_or_empty(entry.get("author")) or "Unknown author"
  modified_text = _annotation_text_or_empty(entry.get("modified_utc")) or "Unknown time"
  preview_lines = _annotation_text_or_empty(entry.get("note")).splitlines()
  preview_text = preview_lines[0].strip() if preview_lines else ""
  if len(preview_text) > max_length:
    preview_text = preview_text[:max_length - 1].rstrip() + "..."
  if preview_text:
    return f"{entry_index}. {author_text} - {modified_text} - {preview_text}"
  return f"{entry_index}. {author_text} - {modified_text}"


def _update_annotation_entry(annotation_entry, author_text, note_text, refresh_function=None):
  global RESIDUE_ANNOTATION_LAST_AUTHOR
  note_fields = _annotation_prepare_note_fields(author_text, note_text)
  if note_fields is None:
    return 0
  normalized_author, normalized_note = note_fields
  RESIDUE_ANNOTATION_LAST_AUTHOR = normalized_author
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
  generic_author_note_entry(
    "Edit residue annotation",
    group["dialog_label"],
    _annotation_text_or_empty(annotation_entry.get("author")) or _annotation_default_author(),
    _annotation_text_or_empty(annotation_entry.get("note")),
    "Save changes",
    lambda author_text, note_text: _update_annotation_entry(
      annotation_entry,
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
  chain_id, resno, ins_code, residue_name_here = group["key"]
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


def _annotation_browser_title(mol_id):
  source_display_name = _annotation_source_display_name(mol_id)
  if source_display_name:
    return (
      f"Residue annotations for molecule #{mol_id}\n"
      f"Source: {source_display_name}; use the widget to export an annotated mmCIF"
    )
  return f"Residue annotations for molecule #{mol_id}\n(No accessible source file path)"


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
  detail_view = Gtk.TextView()
  detail_buffer = detail_view.get_buffer()
  add_button = Gtk.Button(label="Add note for active residue")
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
  detail_view.set_editable(False)
  detail_view.set_cursor_visible(False)
  detail_view.set_wrap_mode(Gtk.WrapMode.WORD_CHAR)
  detail_scrolled.set_child(detail_view)
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
          self.active_residue_timer_id = glib_module.timeout_add(300, self.poll_active_residue_header)
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
      child = list_box.get_first_child()
      while child is not None:
        next_child = child.get_next_sibling()
        list_box.remove(child)
        child = next_child

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
      self.current_index = None
      self.selected_group_key = None
      detail_header.set_text("Selected notes: none")
      detail_buffer.set_text("")
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
      detail_buffer.set_text(_annotation_group_detail_text(group))
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


def _generate_smart_local_extra_restraints_for_mol(mol_id, show_start_message=True):
  if mol_id not in model_molecule_list():
    add_status_bar_text("No active model")
    return None
  if show_start_message:
    add_status_bar_text("Generating smart local restraints...")
    graphics_draw()
  distance_cutoff = 3.7
  max_sequence_separation = 10
  restraint_esd = 0.05
  backbone_donor_names = {"N"}
  backbone_acceptor_names = {"O", "OXT"}

  sequence_index_by_residue = {}
  next_sequence_index_by_chain = {}
  ordered_residue_specs_by_chain = {}
  for residue_spec_list in _all_residue_specs_for_colouring(mol_id):
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
    chain_id, resno, ins_code = residue_spec
    return residue_info_py(mol_id, chain_id, resno, ins_code) or []

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
    is_polymer = _residue_is_polymer(mol_id, residue_spec[0], residue_spec[1], residue_spec[2])
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
    residue_pair_key = ordered_pair_key(record_1["spec"], record_2["spec"])
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
        restraint_key = (
          residue_pair_key,
          (atom_name_1, alt_conf_1),
          (atom_name_2, alt_conf_2)
        )
        reverse_restraint_key = (
          residue_pair_key,
          (atom_name_2, alt_conf_2),
          (atom_name_1, alt_conf_1)
        )
        if restraint_key in seen_restraints or reverse_restraint_key in seen_restraints:
          continue
        distance = atom_distance(atom_1, atom_2)
        if distance is None or distance >= distance_cutoff:
          continue
        add_extra_geman_mcclure_restraint(
          mol_id,
          record_1["spec"][0], record_1["spec"][1], record_1["spec"][2], atom_name_1, alt_conf_1,
          record_2["spec"][0], record_2["spec"][1], record_2["spec"][2], atom_name_2, alt_conf_2,
          distance, restraint_esd
        )
        seen_restraints.add(restraint_key)
        added_restraints += 1

  delete_all_extra_restraints(mol_id)

  added_restraints = 0
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

  all_residue_specs = list(residue_records.keys())
  for residue_spec_1 in all_residue_specs:
    record_1 = residue_records.get(residue_spec_1)
    if not record_1 or not record_1["polymer"]:
      continue
    nearby_residue_specs = residues_near_residue(mol_id, list(residue_spec_1), distance_cutoff) or []
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
  if coot_fitting and hasattr(coot_fitting, "pepflip_active_residue"):
    return coot_fitting.pepflip_active_residue()
  active_atom = closest_atom_simple()
  if not active_atom:
    add_status_bar_text("No active residue")
    return None
  imol = active_atom[0]
  atom_spec = closest_atom_raw()
  chain_id = atom_spec[1]
  res_no = atom_spec[2]
  ins_code = atom_spec[3]
  atom_name = atom_spec[4]
  alt_conf = atom_spec[5]
  if atom_name == " N  ":
    res_no = res_no - 1
  return pepflip(imol, chain_id, res_no, ins_code, alt_conf)


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


def _set_map_radius_both(radius):
  """Keep the standard and EM map-radius settings in sync."""
  set_map_radius(radius)
  set_map_radius_em(radius)


def increase_map_radius():
  current_radius = get_map_radius()
  _set_map_radius_both(current_radius + 2.0)


def decrease_map_radius():
  current_radius = get_map_radius()
  new_radius = max(2.0, current_radius - 2.0)
  _set_map_radius_both(new_radius)


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
  set_draw_map_standard_lines(map_id, 1)
  set_draw_solid_density_surface(map_id, 0)
  if not map_is_difference_map(map_id):
    set_map_colour(
      map_id,
      EM_REFINED_MAP_COLOUR[0],
      EM_REFINED_MAP_COLOUR[1],
      EM_REFINED_MAP_COLOUR[2],
    )
  set_map_material_specular(map_id, 0.5, 64.0)
  set_map_fresnel_settings(map_id, 1, 0.0, 0.3, 3.0)
  set_solid_density_surface_opacity(map_id, 1.0)
  MAP_LOCAL_APPEARANCE_STATE[map_id] = {
    "map_colour": list(map_colour_components_py(map_id)),
    "specular_strength": 0.5,
    "shininess": 64.0,
    "fresnel": (1, 0.0, 0.3, 3.0),
    "opacity": 1.0,
  }


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
    "map_radius": get_map_radius(),
    "clipping_front": get_clipping_plane_front(),
    "clipping_back": get_clipping_plane_back(),
    "map_colour": list(map_colour_components_py(map_id)),
    "opacity": get_solid_density_surface_opacity(map_id),
  }


def _restore_global_view_settings(map_id, fallback_radius):
  saved = MAP_GLOBAL_VIEW_SETTINGS.get(map_id)
  if not saved:
    _set_map_radius_both(fallback_radius)
    return
  _set_map_radius_both(saved["map_radius"])
  set_clipping_front(saved["clipping_front"])
  set_clipping_back(saved["clipping_back"])


def _set_global_view_extent(map_id):
  current_radius = max(get_map_radius(), 1.0)
  global_radius = max(_map_global_extent_radius(map_id), current_radius)
  front = get_clipping_plane_front()
  back = get_clipping_plane_back()
  scale = global_radius / current_radius
  _set_map_radius_both(global_radius)
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
  chain_context = _active_chain_context_or_status()
  if not chain_context:
    return None
  mol_id = chain_context["mol_id"]
  segments=segment_list(mol_id)
  res_here = chain_context["resno"]
  ch_id = chain_context["chain_id"]
  colour_list=[]
  blank_list=[]
  segment_colour=34
  blank_colour=0
  for seg in segments:
    if (res_here>=seg[2]) and (res_here<=seg[3]) and (ch_id==seg[1]):
      res_start=seg[2]
      res_end=seg[3]
      ch_id=seg[1]
      for res in range(res_start,res_end+1):
        colour_list.append(([ch_id,res,""],segment_colour))
    else:
      res_start=seg[2]
      res_end=seg[3]
      ch_id_here=seg[1]
      for res in range(res_start,res_end+1):
        blank_list.append(([ch_id_here,res,""],blank_colour))
  _apply_user_defined_residue_colours(mol_id, blank_list, colour_list)


def _all_residue_specs_for_colouring(mol_id):
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


def _clear_custom_colour_additional_representations(mol_id):
  handles = CUSTOM_COLOUR_ADDITIONAL_REPRESENTATIONS.pop(mol_id, [])
  for handle in handles:
    try:
      delete_additional_representation(mol_id, handle)
    except Exception:
      pass
  try:
    clear_user_defined_atom_colours(mol_id)
  except Exception:
    pass


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
  return coot.set_user_defined_atom_colour_py(mol_id, atom_colour_list)


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
  _clear_custom_colour_additional_representations(mol_id)
  add_status_bar_text("Cleared custom colour representations")


def _apply_user_defined_residue_colours(mol_id, blank_res_list, colour_list, info_message=None):
  del blank_res_list
  ensure_user_defined_colour_table()
  _clear_custom_colour_additional_representations(mol_id)
  if not colour_list:
    if info_message:
      info_dialog(info_message)
    return None
  clear_user_defined_atom_colours(mol_id)
  colour_count = _set_user_defined_atom_colour_by_residue_atoms_py(mol_id, colour_list)
  if not colour_count:
    _clear_custom_colour_additional_representations(mol_id)
    info_dialog("Failed to set user-defined atom colours for the requested residues.")
    return None
  handles = []
  seen_residue_specs = set()
  for residue_spec, _colour_index in colour_list:
    residue_key = tuple(residue_spec[:3])
    if residue_key in seen_residue_specs:
      continue
    seen_residue_specs.add(residue_key)
    handle = _make_custom_colour_additional_representation(mol_id, residue_spec)
    if handle is not None:
      handles.append(handle)
  CUSTOM_COLOUR_ADDITIONAL_REPRESENTATIONS[mol_id] = handles
  clear_user_defined_atom_colours(mol_id)
  if not handles:
    info_dialog("Failed to create a custom-colour additional representation.")
    return None
  if info_message:
    info_dialog(info_message)


def color_by_rama_native(mol_id):
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
  _apply_user_defined_residue_colours(
    mol_id,
    blank_res_list,
    rama_colour_list,
    "Ramachandran coloring:\n\nRed = outlier (<0.2%)\n\nOrange = allowed/disfavored (<2%)",
  )


def color_by_rama_native_for_active_residue():
  mol_id=_active_polymer_molecule_for_colouring("Ramachandran coloring")
  if mol_id is None:
    return None
  return color_by_rama_native(mol_id)


def color_by_density_fit_native(mol_id):
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
    if not isinstance(item, list) or len(item) < 2:
      continue
    residue_spec = _coot_residue_spec_from_spec(item[0])
    score=item[1]
    if residue_spec is None:
      continue
    if score < 0.0:
      score=0.0
    if score > 1.0:
      score=1.0
    colour_index=int((1.0-score)*31+2)
    density_colour_list.append((list(residue_spec), colour_index))
  _apply_user_defined_residue_colours(
    mol_id,
    blank_res_list,
    density_colour_list,
    "Active molecule colored by model/map correlation, in spectral coloring (blue=CC 1.0, red=CC 0.0)",
  )


def color_by_density_fit_native_for_active_residue():
  mol_id=_active_polymer_molecule_for_colouring("Density-fit coloring")
  if mol_id is None:
    return None
  return color_by_density_fit_native(mol_id)


def color_by_ncs_difference(mol_id):
  if mol_id not in model_molecule_list():
    info_dialog("You need an active model for NCS-difference coloring.")
    return None
  ncs_data=None
  for chain_id in chain_ids(mol_id):
    try:
      diffs=ncs_chain_differences(mol_id, chain_id)
    except Exception:
      diffs=False
    if diffs:
      ncs_data=diffs
      break
  if not ncs_data:
    info_dialog("No NCS-difference data were found for the active molecule.")
    return None

  blank_colour=0
  blank_res_list=[]
  ncs_scores={}

  for residue_spec in _all_residue_specs_for_colouring(mol_id):
    blank_res_list.append((residue_spec, blank_colour))

  for i in range(0, len(ncs_data), 3):
    try:
      peer_chain_id=ncs_data[i]
      current_target_chain_id=ncs_data[i+1]
      residue_diffs=ncs_data[i+2]
    except Exception:
      continue
    if not isinstance(residue_diffs, list):
      continue
    for residue_diff in residue_diffs:
      if not isinstance(residue_diff, list) or len(residue_diff) < 3:
        continue
      peer_residue=residue_diff[0]
      target_residue=residue_diff[1]
      mean_diff=residue_diff[2]
      try:
        mean_diff=float(mean_diff)
      except Exception:
        continue
      if isinstance(peer_residue, list) and len(peer_residue) >= 2:
        peer_spec=(peer_chain_id, peer_residue[0], peer_residue[1])
        ncs_scores[peer_spec]=max(ncs_scores.get(peer_spec, 0.0), mean_diff)
      if isinstance(target_residue, list) and len(target_residue) >= 2:
        target_spec=(current_target_chain_id, target_residue[0], target_residue[1])
        ncs_scores[target_spec]=max(ncs_scores.get(target_spec, 0.0), mean_diff)

  if not ncs_scores:
    info_dialog("No NCS-difference values were available for coloring.")
    return None

  ncs_colour_list=[]
  for residue_spec, mean_diff in ncs_scores.items():
    normalized_score=mean_diff/2.0
    if normalized_score < 0.0:
      normalized_score=0.0
    if normalized_score > 1.0:
      normalized_score=1.0
    colour_index=int(normalized_score*31+2)
    ncs_colour_list.append(([residue_spec[0], residue_spec[1], residue_spec[2]], colour_index))

  _apply_user_defined_residue_colours(
    mol_id,
    blank_res_list,
    ncs_colour_list,
    "Active molecule colored by NCS difference, in spectral coloring (blue=low difference, red=high difference)",
  )


def color_by_ncs_difference_for_active_residue():
  mol_id=_active_polymer_molecule_for_colouring("NCS-difference coloring")
  if mol_id is None:
    return None
  return color_by_ncs_difference(mol_id)


def color_by_clash_score(mol_id):
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
    colour_index=int(normalized_score*31+2)
    clash_colour_list.append(([residue_spec[0], residue_spec[1], residue_spec[2]], colour_index))

  _apply_user_defined_residue_colours(
    mol_id,
    blank_res_list,
    clash_colour_list,
    "Active molecule colored by per-residue maximum clash overlap, in spectral coloring (blue=low clash, red=high clash)",
  )


def color_by_clash_score_for_active_molecule():
  mol_id=_active_molecule_or_status()
  if mol_id is None:
    return None
  return color_by_clash_score(mol_id)

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
REPRESENTATION_SEQUENCE = [
  graphics_to_ca_plus_ligands_representation,
  graphics_to_ca_plus_ligands_and_sidechains_representation,
  graphics_to_ca_plus_ligands_sec_struct_representation,
  graphics_to_rainbow_representation,
  graphics_to_bonds_representation,
  graphics_to_b_factor_representation,
]
DEFAULT_REPRESENTATION_INDEX = 4
cycle_rep_flag = {}

def _cycle_rep_flag_or_default(mol_id):
  flag = cycle_rep_flag.get(mol_id, DEFAULT_REPRESENTATION_INDEX)
  if flag not in range(len(REPRESENTATION_SEQUENCE)):
    flag = DEFAULT_REPRESENTATION_INDEX
  cycle_rep_flag[mol_id] = flag
  return flag

def _apply_representation_index(mol_id, flag):
  REPRESENTATION_SEQUENCE[flag](mol_id)
  cycle_rep_flag[mol_id] = flag
  return flag

def cycle_rep_up(mol_id,flag=None):
  current_flag = _cycle_rep_flag_or_default(mol_id) if flag is None else flag
  next_flag = (current_flag + 1) % len(REPRESENTATION_SEQUENCE)
  return _apply_representation_index(mol_id, next_flag)


    
def cycle_rep_down(mol_id,flag=None):
  current_flag = _cycle_rep_flag_or_default(mol_id) if flag is None else flag
  next_flag = (current_flag - 1) % len(REPRESENTATION_SEQUENCE)
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
  
def redo_visible():
  residue = _active_residue_or_status()
  if not residue:
    return None
  set_undo_molecule(residue[0])
  apply_redo()

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
  turn_off_backup(mol_id)
  try:
    update_go_to_atom_from_current_position()
    if rotamer_number>=n_rots:
      rotamer_number=0
      set_residue_to_rotamer_number(mol_id,ch_id,res_here,ins_code,alt_conf,rotamer_number)
    else:
      rotamer_number=rotamer_number+1
      set_residue_to_rotamer_number(mol_id,ch_id,res_here,ins_code,alt_conf,rotamer_number)
  finally:
    turn_on_backup(mol_id)
  
#Toggle environment distances
def toggle_env_dist():
  if show_environment_distances_state()==1:
    set_show_environment_distances(0)
  else:
    set_show_environment_distances(1)
    
#Check if residue is in a polymer by comparing resname to list
def is_polymer_residue(mol_id,ch_id,sn):
  valid_resnames=['A','C','T','G','U','ALA','UNK','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','ILE','LEU','LYS','MET','MSE','PHE','PRO','SER','THR','TRP','TYR','VAL']
  resname=resname_from_serial_number(mol_id,ch_id,sn)
  if resname in valid_resnames:
    return 1
  else:
    return 0
    
#Return last residue in polymer
def last_polymer_residue(mol_id,ch_id):
  if ((valid_model_molecule_qm(mol_id)) and (ch_id in chain_ids(mol_id))):
    if not is_solvent_chain_qm(mol_id,ch_id):
      n=chain_n_residues(ch_id,mol_id)-1
      valid_resnames=['A','C','T','G','U','ALA','UNK','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','ILE','LEU','LYS','MET','MSE','PHE','PRO','SER','THR','TRP','TYR','VAL']
      while resname_from_serial_number(mol_id,"%s"%(ch_id),n) not in valid_resnames:
        n=n-1 
      result=seqnum_from_serial_number(mol_id,"%s"%(ch_id),n)
      return result
  else:
    return -1
    
def first_polymer_residue(mol_id,ch_id):
  if ((valid_model_molecule_qm(mol_id)) and (ch_id in chain_ids(mol_id))):
    if not is_solvent_chain_qm(mol_id,ch_id):
      n=0
      valid_resnames=['A','C','T','G','U','ALA','UNK','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','ILE','LEU','LYS','MET','MSE','PHE','PRO','SER','THR','TRP','TYR','VAL']
      while resname_from_serial_number(mol_id,"%s"%(ch_id),n) not in valid_resnames:
        n=n+1 
      result=seqnum_from_serial_number(mol_id,"%s"%(ch_id),n)
      return result
  else:
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
def is_last_polymer_residue_sn(mol_id,ch_id,sn):
  if ((valid_model_molecule_qm(mol_id)) and (ch_id in chain_ids(mol_id))):
    if not is_solvent_chain_qm(mol_id,ch_id):
      valid_resnames=['A','C','T','G','U','ALA','UNK','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','ILE','LEU','LYS','MET','MSE','PHE','PRO','SER','THR','TRP','TYR','VAL']
      if (resname_from_serial_number(mol_id,"%s"%(ch_id),sn) in valid_resnames) and (resname_from_serial_number(mol_id,"%s"%(ch_id),sn+1) not in valid_resnames):
        return 1
  else:
    return -1

#get serial number from resnum
def get_sn_from_resno(mol_id,ch_id,resno):
  sn=chain_n_residues(ch_id,mol_id)-1
  sn2=0
  resno_out=""
  resno_out_2=""
  if resno>last_residue(mol_id,ch_id) or resno<first_residue(mol_id,ch_id):
    return -1
  elif does_residue_exist_p(mol_id,ch_id,resno,"")==0:
    return -1
  elif resno==first_residue(mol_id,ch_id):
    return 0
  elif resno==chain_n_residues(ch_id,mol_id):
    return chain_n_residues(ch_id,mol_id)
  else:
    while (resno_out!=resno) and (resno_out_2!=resno):
      resno_out=seqnum_from_serial_number(mol_id,ch_id,sn)
      resno_out_2=seqnum_from_serial_number(mol_id,ch_id,sn2)
      sn=sn-1
      sn2=sn+1
    if resno_out_2==resno:
      return int(sn2+1)
    else:
      return int(sn+1)

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
def is_term_type_mc(mol_id,ch_id,resno):
  sn=get_sn_from_resno(mol_id,ch_id,resno)
  if (type(sn) is int) and (is_polymer_residue(mol_id,ch_id,sn)==1) and (is_polymer_residue(mol_id,ch_id,sn+1)==1):
    resn_here=resno
    resn_next=seqnum_from_serial_number(mol_id,ch_id,sn+1)
    diff_next=resn_next-resn_here
    if (diff_next>=2):
      return 1
    else:
      return 0
  else:
    return 0
    

#check if res is at C-term side of break in mid chain
def is_term_type_mc_sn(mol_id,ch_id,sn):
  if sn in range(1,chain_n_residues(ch_id,mol_id)):
    resn_here=seqnum_from_serial_number(mol_id,ch_id,sn)
    if (is_polymer_residue(mol_id,ch_id,sn)==1) and (is_polymer_residue(mol_id,ch_id,sn+1)==1):
      resn_next=seqnum_from_serial_number(mol_id,ch_id,sn+1)
      diff_next=resn_next-resn_here
      if (diff_next>=2):
        return 1
      else:
        return 0
    else:
      return 0
  else:
    return 0
    
#check if res is at N-term side of break in mid chain
def is_term_type_mn(mol_id,ch_id,resno):
  sn=get_sn_from_resno(mol_id,ch_id,resno)
  if type(sn) is int:
    if sn in range(1,chain_n_residues(ch_id,mol_id)) and (is_polymer_residue(mol_id,ch_id,sn)==1):
      resn_here=resno
      resn_prev=seqnum_from_serial_number(mol_id,ch_id,sn-1)
      diff_prev=resn_here-resn_prev
      if (diff_prev>=2):
        return 1
      else:
        return 0
    else:
      return 0
  else:
    return 0
    
def is_term_type_mn_sn(mol_id,ch_id,sn):
  if sn in range(1,chain_n_residues(ch_id,mol_id)):
    resn_here=seqnum_from_serial_number(mol_id,ch_id,sn) 
    if (is_polymer_residue(mol_id,ch_id,sn)==1):
      resn_prev=seqnum_from_serial_number(mol_id,ch_id,sn-1)
      diff_prev=resn_here-resn_prev
      if (diff_prev>=2):
        return 1
      else:
        return 0
    else:
      return 0
  else:
    return 0
  
        
def first_residue_in_seg(mol_id,ch_id,resno):
  sn_here=get_sn_from_resno(mol_id,ch_id,resno)
  sn=sn_here
  while is_term_type_mn_sn(mol_id,ch_id,sn)==0 and sn>0 and seqnum_from_serial_number(mol_id,ch_id,sn-1)!=-10000 and is_polymer_residue(mol_id,ch_id,sn-1)==1:
    sn=sn-1
  res_start=seqnum_from_serial_number(mol_id,ch_id,sn)
  return res_start
  
def last_residue_in_seg(mol_id,ch_id,resno):
  sn_here=get_sn_from_resno(mol_id,ch_id,resno)
  sn=sn_here
  while is_term_type_mc_sn(mol_id,ch_id,sn)==0 and seqnum_from_serial_number(mol_id,ch_id,sn+1)!=-10000 and is_polymer_residue(mol_id,ch_id,sn+1)==1:
    sn=sn+1
  res_end=seqnum_from_serial_number(mol_id,ch_id,sn)
  return res_end


#Get serial number of active residue
def sn_of_active_res():
  chain_context = _active_chain_context_or_status()
  if not chain_context:
    return None
  mol_id = chain_context["mol_id"]
  ch_id = chain_context["chain_id"]
  resn_to_match = chain_context["resno"]
  for sn in range(chain_n_residues(ch_id, mol_id) + 1):
    if seqnum_from_serial_number(mol_id, "%s" % ch_id, sn) == resn_to_match:
      return sn
  return None
    
#Get monomer and delete hydrogens
def get_monomer_no_H(mon):
  get_monomer(mon)
  delete_hydrogens(molecule_number_list()[-1])    
    
#Return list of segments in active mol
def segment_list(mol_id):
  sn=0
  list_out=[]
  for ch_id in chain_ids(mol_id):
    while is_polymer_residue(mol_id,ch_id,sn+1)==1:
      if sn==0:
        res_start=seqnum_from_serial_number(mol_id,ch_id,sn)
        while is_term_type_mc_sn(mol_id,ch_id,sn)==0 and seqnum_from_serial_number(mol_id,ch_id,sn+1)!=-10000 and is_polymer_residue(mol_id,ch_id,sn+1)==1:
          sn=sn+1
        res_end=seqnum_from_serial_number(mol_id,ch_id,sn)
        list_out.append([mol_id,ch_id,res_start,res_end])
      else:
        while is_term_type_mn_sn(mol_id,ch_id,sn)==0:
          sn=sn+1
        res_start=seqnum_from_serial_number(mol_id,ch_id,sn)
        while is_term_type_mc_sn(mol_id,ch_id,sn)==0 and seqnum_from_serial_number(mol_id,ch_id,sn+1)!=-10000 and is_polymer_residue(mol_id,ch_id,sn+1)==1:
          sn=sn+1
        res_end=seqnum_from_serial_number(mol_id,ch_id,sn)
        list_out.append([mol_id,ch_id,res_start,res_end])
    sn=0
  return list_out
  
def segment_list_chain(mol_id,ch_id):
  sn=0
  list_out=[]
  while is_polymer_residue(mol_id,ch_id,sn+1)==1:
    if sn==0:
      res_start=seqnum_from_serial_number(mol_id,ch_id,sn)
      while is_term_type_mc_sn(mol_id,ch_id,sn)==0 and seqnum_from_serial_number(mol_id,ch_id,sn+1)!=-10000 and is_polymer_residue(mol_id,ch_id,sn+1)==1:
        sn=sn+1
      res_end=seqnum_from_serial_number(mol_id,ch_id,sn)
      list_out.append([mol_id,ch_id,res_start,res_end])
    else:
      while is_term_type_mn_sn(mol_id,ch_id,sn)==0:
        sn=sn+1
      res_start=seqnum_from_serial_number(mol_id,ch_id,sn)
      while is_term_type_mc_sn(mol_id,ch_id,sn)==0 and seqnum_from_serial_number(mol_id,ch_id,sn+1)!=-10000 and is_polymer_residue(mol_id,ch_id,sn+1)==1:
        sn=sn+1
      res_end=seqnum_from_serial_number(mol_id,ch_id,sn)
      list_out.append([mol_id,ch_id,res_start,res_end])
  return list_out
  
    
#Local cylinder refinement around the active residue:
# main polymer range +/- 5 residues, contact-expanded secondary ranges,
# then local polymer gap-filling/pruning while keeping non-polymer contacts.
def auto_refine():
  if imol_refinement_map()==-1:
    info_dialog("You must set a refinement map!")
    return None
  residue = _active_residue_or_status()
  if not residue:
    return None
  mol_id=residue[0]
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

  if _residue_is_polymer(mol_id, ch_id, resno, active_ins_code):
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
  previous_immediate_replacement = refinement_immediate_replacement_state()
  try:
    set_refinement_immediate_replacement(0)
    refine_residues_py(mol_id, residue_list)
    add_status_bar_text("Local cylinder refinement started")
  finally:
    set_refinement_immediate_replacement(previous_immediate_replacement)
   

#**** "Custom menu item functions ****
#Fits all polymer chains to map
def rigid_fit_all_chains():
  mol_id = _active_molecule_or_status()
  if mol_id is None:
    return None
  turn_off_backup(mol_id)
  try:
    for ch_id in chain_ids(mol_id):
      if is_polymer(mol_id, ch_id) == 1:
        rigid_body_refine_by_atom_selection(mol_id, "//%s//" % ch_id)
        accept_regularizement()
  finally:
    turn_on_backup(mol_id)
  
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
  new_molecule_by_atom_selection(mol_id, "//%s//" % ch_id)
  turn_off_backup(mol_id)
  try:
    while (is_polymer(mol_id, ch_id) == 1) or (is_solvent_chain_p(mol_id, ch_id) != -1):
      first_res = first_residue(mol_id, ch_id)
      last_res = last_residue(mol_id, ch_id)
      delete_residue_range(mol_id, ch_id, first_res, last_res)
  finally:
    turn_on_backup(mol_id)

#Faster rigid body fit (not relevant for smaller ranges, but much faster for larger ranges)
def fast_rigid_fit(res_start,res_end,ch_id,mol_id):
  if (mol_id in model_molecule_number_list()) and (ch_id in chain_ids(mol_id)):
    turn_off_backup(mol_id)
    try:
      ins_code=""
      sn_max=chain_n_residues(ch_id,mol_id)-1
      res_min=seqnum_from_serial_number(mol_id,ch_id,0)
      res_max=seqnum_from_serial_number(mol_id,ch_id,sn_max)
      if res_start>res_end:
        res_start,res_end=res_end,res_start
      while not residue_exists_qm(mol_id,ch_id,res_start,ins_code) and res_start<=res_max:
       res_start=res_start+1
      while not residue_exists_qm(mol_id,ch_id,res_end,ins_code) and res_end>=res_min:
       res_end=res_end-1
      if res_start<=res_end:
        new_molecule_by_atom_selection(mol_id, "//{ch_id}/{res_start}-{res_end}/".format(ch_id=ch_id,res_start=res_start,res_end=res_end))
        mol_id_new=model_molecule_number_list()[-1]
        rigid_body_refine_by_atom_selection(mol_id_new,"/ /")
        accept_regularizement()
        delete_residue_range(mol_id,ch_id,res_start,res_end) #delete copied range from original mol
        merge_molecules([mol_id_new],mol_id) #Merge fit segment back into original mol
        change_chain_id_with_result(mol_id,chain_ids(mol_id)[-1],ch_id,1,res_start,res_end) #Merge chains
        close_molecule(mol_id_new)
    finally:
      turn_on_backup(mol_id)
  else:
    info_dialog("The specified chain or molecule does not exist!")

#Copy active segment
def copy_active_segment():
  segment_context = _active_segment_context_or_status()
  if not segment_context:
    return None
  new_molecule_by_atom_selection(
    segment_context["mol_id"],
    "//{ch_id}/{res_start}-{res_end}/".format(
      ch_id=segment_context["chain_id"],
      res_start=segment_context["segment_start"],
      res_end=segment_context["segment_end"],
    ),
  )

#Cut active segment
def cut_active_segment():
  segment_context = _active_segment_context_or_status()
  if not segment_context:
    return None
  mol_id = segment_context["mol_id"]
  ch_id = segment_context["chain_id"]
  res_start = segment_context["segment_start"]
  res_end = segment_context["segment_end"]
  new_molecule_by_atom_selection(mol_id, "//{ch_id}/{res_start}-{res_end}/".format(ch_id=ch_id, res_start=res_start, res_end=res_end))
  delete_residue_range(mol_id, ch_id, res_start, res_end)

#Delete active segment
def delete_active_segment():
  segment_context = _active_segment_context_or_status()
  if not segment_context:
    return None
  delete_residue_range(
    segment_context["mol_id"],
    segment_context["chain_id"],
    segment_context["segment_start"],
    segment_context["segment_end"],
  )


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
  sharpen(imol_refinement_map(),200)
  try:
    fit_chain_to_map_by_random_jiggle(chain_context["mol_id"], chain_context["chain_id"], 1000, 0.1)
  finally:
    sharpen(imol_refinement_map(),0)
    
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
  turn_off_backup(mol_id)
  try:
    for ch_id in chain_ids(mol_id):
      if is_polymer(mol_id,ch_id)==1:
        fit_chain_to_map_by_random_jiggle(mol_id,ch_id,1000,0.1)
        accept_regularizement()
  finally:
    turn_on_backup(mol_id)

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
  user_defined_click(2,fit_polyala)
  
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
    if resno_2>=resno_1:
      refine_zone(mol_id_1,ch_id_1,resno_1,resno_2,alt_conf_1)
    else:
      refine_zone(mol_id_1,ch_id_1,resno_2,resno_1,alt_conf_1)
  user_defined_click(2,refine_click_a)
  
#Refine a clicked residue range.
def refine_residues_click():
  add_status_bar_text("Click 2 atoms to refine zone")
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
  user_defined_click(2,refine_residues_click_a)
  
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
  if (mol_id_1!=mol_id_2) or (ch_id_1!=ch_id_2):
    info_dialog(error_message)
    return None
  if resno_1 == resno_2 and not allow_single_residue:
    info_dialog(error_message)
    return None
  return (mol_id_1, ch_id_1, min(resno_1, resno_2), max(resno_1, resno_2))


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
    atom_sel = _atom_selection_for_residue_range(ch_id, res_start, res_end)
    new_molecule_by_atom_selection(mol_id, atom_sel)
  user_defined_click(2,copy_frag)
  
#Cut fragment (click start and end)
def cut_frag_by_click():
  def cut_frag(res1,res2):
    clicked_range = _clicked_same_chain_residue_range(res1, res2)
    if clicked_range is None:
      return None
    mol_id, ch_id, res_start, res_end = clicked_range
    atom_sel = _atom_selection_for_residue_range(ch_id, res_start, res_end)
    turn_off_backup(mol_id)
    try:
      new_molecule_by_atom_selection(mol_id, atom_sel)
      if res_start == res_end:
        delete_residue_range(mol_id, ch_id, res_start)
      else:
        delete_residue_range(mol_id, ch_id, res_start, res_end)
    finally:
      turn_on_backup(mol_id)
  user_defined_click(2,cut_frag)


#Delete sidechain range (click start and end)
def delete_sidechain_range(mol_id, ch_id, res_start, res_end):
  for resno in range(res_start, res_end + 1):
    coot.delete_residue_sidechain(mol_id, ch_id, resno, "", 0)

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
    turn_off_backup(mol_id)
    try:
      mutate_residue_range(mol_id,ch_id,res_start,res_end,target_seq)
      for resno in range(res_start,res_end+1):
        set_residue_name(mol_id,ch_id,resno,"","UNK")
    finally:
      turn_on_backup(mol_id)
  user_defined_click(2,mutate_residue_range_by_click_b)

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
    turn_off_backup(mol_id)
    try:
      mutate_residue_range(mol_id,ch_id,res_start,res_end,target_seq)
      for resno in range(res_start,res_end+1):
        set_residue_name(mol_id,ch_id,resno,"","ALA")
    finally:
      turn_on_backup(mol_id)
  user_defined_click(2,mutate_residue_range_by_click_ala_b)
  
#Force addition of residue - useful when
#Coot says "No acceptable position found"
# but density is clear.
def force_add_terminal_residue():
  def force_addition(res1):
    mol_id=_click_spec_imol(res1)
    ch_id=_click_spec_chain_id(res1)
    res_no=_click_spec_res_no(res1)
    ins_code=_click_spec_ins_code(res1)
    res_type="auto"
    add_terminal_residue_using_phi_psi(mol_id,ch_id,res_no,
    res_type,-57.82,-47)
    if residue_exists_qm(mol_id,ch_id,res_no+1,ins_code):
      set_b_factor_residue_range(mol_id,ch_id,res_no+1,res_no+1,default_new_atoms_b_factor())
    elif residue_exists_qm(mol_id,ch_id,res_no-1,ins_code):
      set_b_factor_residue_range(mol_id,ch_id,res_no-1,res_no-1,default_new_atoms_b_factor())
    sort_residues(mol_id)
  user_defined_click(1,force_addition)
  
def force_add_terminal_residue_noclick(mol_id,ch_id,res_no):
  res_type="auto"
  add_terminal_residue_using_phi_psi(mol_id,ch_id,res_no,
  res_type,-57.82,-47)
  if residue_exists_qm(mol_id,ch_id,res_no+1,""):
    set_b_factor_residue_range(mol_id,ch_id,res_no+1,res_no+1,default_new_atoms_b_factor())
  elif residue_exists_qm(mol_id,ch_id,res_no-1,""):
    set_b_factor_residue_range(mol_id,ch_id,res_no-1,res_no-1,default_new_atoms_b_factor())
  sort_residues(mol_id)

def force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,res_no,phi,psi):
  res_type="auto"
  add_terminal_residue_using_phi_psi(mol_id,ch_id,res_no,
  res_type,float(phi),float(psi))
  if residue_exists_qm(mol_id,ch_id,res_no+1,""):
    set_b_factor_residue_range(mol_id,ch_id,res_no+1,res_no+1,default_new_atoms_b_factor())
  elif residue_exists_qm(mol_id,ch_id,res_no-1,""):
    set_b_factor_residue_range(mol_id,ch_id,res_no-1,res_no-1,default_new_atoms_b_factor())
  sort_residues(mol_id)

def force_add_terminal_residue_noclick_strand(mol_id,ch_id,res_no):
  res_type="auto"
  add_terminal_residue_using_phi_psi(mol_id,ch_id,res_no,
  res_type,-139,135)
  if residue_exists_qm(mol_id,ch_id,res_no+1,""):
    set_b_factor_residue_range(mol_id,ch_id,res_no+1,res_no+1,default_new_atoms_b_factor())
  elif residue_exists_qm(mol_id,ch_id,res_no-1,""):
    set_b_factor_residue_range(mol_id,ch_id,res_no-1,res_no-1,default_new_atoms_b_factor())
  sort_residues(mol_id)

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
  next_cycle_index = (cycle_index + 1) % len(RESIDUE_TORSION_CYCLE_VALUES)
  if vary_phi:
    phi = chosen_angle
    psi = current_psi
    residue_phi_cycle = next_cycle_index
  else:
    phi = current_phi
    psi = chosen_angle
    residue_psi_cycle = next_cycle_index
  set_go_to_atom_molecule(mol_id)
  set_new_atom_b_fac_to_mean()
  set_go_to_atom_chain_residue_atom_name(ch_id, anchor_resno, atom_name)
  delete_residue(mol_id, ch_id, anchor_resno, "")
  force_add_terminal_residue_noclick_phi_psi(mol_id, ch_id, insertion_resno, phi, psi)
  sort_residues(mol_id)
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
  turn_off_backup(mol_id)
  try:
    for ch_id in chain_ids(mol_id):
      for resn in range(0, chain_n_residues(ch_id, mol_id)):
        if resname_from_serial_number(mol_id, ch_id, resn) == source_resname:
          seqnum=seqnum_from_serial_number(mol_id,ch_id,resn)
          ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,resn))
          delete_residue_with_full_spec(mol_id,1,ch_id,seqnum,ins_id,"B")
          delete_residue_with_full_spec(mol_id,1,ch_id,seqnum,ins_id,"C")
          mutate(mol_id,ch_id,seqnum,ins_id,target_resname)
          if (map_id!=-1):
            auto_fit_best_rotamer(seqnum,"",ins_id,ch_id,mol_id,map_id,1,0.01)
  finally:
    turn_on_backup(mol_id)


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
  renumber_residue_range(mol_id,ch_id,first_res,resn,1)
  _refresh_extra_restraints_display(mol_id)
  r1=resn-1
  r2=resn+2
  set_refinement_immediate_replacement(1)
  refine_zone(mol_id,ch_id,r1,r2,"")
  accept_regularizement()
  set_refinement_immediate_replacement(0)

#Lengthen loop by one residue
def lengthen_loop():
  active_atom=_active_residue_or_status()
  if not active_atom:
    return None
  mol_id=active_atom[0]
  ch_id=active_atom[1]
  resn=active_atom[2]
  first_res=seqnum_from_serial_number(mol_id,"%s"%(ch_id),0)
  renumber_residue_range(mol_id,ch_id,first_res,resn,-1)
  sort_residues(mol_id)
  _refresh_extra_restraints_display(mol_id)
  r1=resn-1
  r2=resn
  set_refinement_immediate_replacement(1)
  try:
    status = add_residue_by_map_fit(mol_id,ch_id,r1,"auto",1)
    if not status:
      info_dialog("Failed to add a residue into the loop gap")
      return None
    refine_zone(mol_id,ch_id,r1,r2+1,"")
    accept_regularizement()
  finally:
    set_refinement_immediate_replacement(0)

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
def color_rotamer_outliers_and_missing_atoms(mol_id):
  missing_atoms_list=[]
  missing_atoms_colour=2
  rotamer_outlier_list=[]
  rotamer_outlier_colour=34
  blank_colour=0
  for x in missing_atom_info(mol_id):
    missing_atoms_spec=[(x,missing_atoms_colour)]
    missing_atoms_list=missing_atoms_list+missing_atoms_spec
  for ch_id in chain_ids(mol_id):
    first_res=first_residue(mol_id,ch_id)
    last_res=last_residue(mol_id,ch_id)
    for resn in range(first_res,last_res):
      if residue_exists_qm(mol_id,ch_id,resn,""):
        rot_prob=rotamer_score(mol_id,ch_id,resn,"","")
        if rot_prob<0.5 and rot_prob>0.0:
          rotamer_outlier_spec=[([ch_id,resn,""],rotamer_outlier_colour)]
          rotamer_outlier_list=rotamer_outlier_list+rotamer_outlier_spec
        else:
          rotamer_outlier_spec=[([ch_id,resn,""],blank_colour)]
          rotamer_outlier_list=rotamer_outlier_list+rotamer_outlier_spec
  colour_list = rotamer_outlier_list + missing_atoms_list
  _apply_user_defined_residue_colours(mol_id, [], colour_list)


def color_polars_and_hphobs(mol_id):
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
  for ch_id in chain_ids(mol_id):
    sn_max=chain_n_residues(ch_id,mol_id)
    for sn in range(0,sn_max+1):
      resname_here=resname_from_serial_number(mol_id,ch_id,sn)
      if resname_here in polar_list:
        resn=seqnum_from_serial_number(mol_id,ch_id,sn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
        residue_to_color=[([ch_id,resn,ins_id],polar_colour)]
        polar_res_list=polar_res_list+residue_to_color
      elif resname_here in hphob_list:
        resn=seqnum_from_serial_number(mol_id,ch_id,sn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
        residue_to_color=[([ch_id,resn,ins_id],hphob_colour)]
        hphob_res_list=hphob_res_list+residue_to_color
      elif resname_here=="GLY":
        resn=seqnum_from_serial_number(mol_id,ch_id,sn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
        residue_to_color=[([ch_id,resn,ins_id],gly_colour)]
        gly_res_list=gly_res_list+residue_to_color
      elif resname_here=="PRO":
        resn=seqnum_from_serial_number(mol_id,ch_id,sn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
        residue_to_color=[([ch_id,resn,ins_id],pro_colour)]
        pro_res_list=pro_res_list+residue_to_color
      else:
        resn=seqnum_from_serial_number(mol_id,ch_id,sn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
        residue_to_color=[([ch_id,resn,ins_id],blank_colour)]
        blank_res_list=blank_res_list+residue_to_color
  colour_list = polar_res_list + hphob_res_list + gly_res_list + pro_res_list
  _apply_user_defined_residue_colours(mol_id, blank_res_list, colour_list)
  
def color_by_charge(mol_id):
  pos_list=["ARG","LYS","HIS"]
  neg_list=["GLU","ASP"]
  pos_colour=4
  neg_colour=31
  blank_colour=0
  pos_res_list=[]
  neg_res_list=[]
  blank_res_list=[]
  for ch_id in chain_ids(mol_id):
    sn_max=chain_n_residues(ch_id,mol_id)
    for sn in range(0,sn_max+1):
      resname_here=resname_from_serial_number(mol_id,ch_id,sn)
      if resname_here in pos_list:
        resn=seqnum_from_serial_number(mol_id,ch_id,sn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
        residue_to_color=[([ch_id,resn,ins_id],pos_colour)]
        pos_res_list=pos_res_list+residue_to_color
      elif resname_here in neg_list:
        resn=seqnum_from_serial_number(mol_id,ch_id,sn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
        residue_to_color=[([ch_id,resn,ins_id],neg_colour)]
        neg_res_list=neg_res_list+residue_to_color
      else:
        resn=seqnum_from_serial_number(mol_id,ch_id,sn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
        residue_to_color=[([ch_id,resn,ins_id],blank_colour)]
        blank_res_list=blank_res_list+residue_to_color
  colour_list = pos_res_list + neg_res_list
  _apply_user_defined_residue_colours(mol_id, blank_res_list, colour_list)


def _append_chain_residue_colours(colour_rows, mol_id, ch_id, colour_index):
  sn_max = chain_n_residues(ch_id, mol_id)
  for sn in range(0, sn_max+1):
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
    
def color_protein_na(mol_id):
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
  _apply_user_defined_residue_colours(mol_id, blank_res_list, colour_list)


def color_waters(mol_id):
  water_colour=31
  blank_colour=0
  water_list=[]
  blank_res_list=[]
  for ch_id in chain_ids(mol_id):
    sn_max=chain_n_residues(ch_id,mol_id)
    for sn in range(0,sn_max+1):
      resname_here=resname_from_serial_number(mol_id,ch_id,sn)
      resn=seqnum_from_serial_number(mol_id,ch_id,sn)
      ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
      colour_entry = ([ch_id,resn,ins_id], water_colour if resname_here=="HOH" else blank_colour)
      if resname_here=="HOH":
        water_list.append(colour_entry)
      else:
        blank_res_list.append(colour_entry)
  _apply_user_defined_residue_colours(mol_id, blank_res_list, water_list)


def color_rotamer_outliers_and_missing_atoms_for_active_molecule():
  mol_id = _active_molecule_or_status()
  if mol_id is None:
    return None
  return color_rotamer_outliers_and_missing_atoms(mol_id)


def color_polars_and_hphobs_for_active_molecule():
  mol_id = _active_molecule_or_status()
  if mol_id is None:
    return None
  return color_polars_and_hphobs(mol_id)


def color_by_charge_for_active_molecule():
  mol_id = _active_molecule_or_status()
  if mol_id is None:
    return None
  return color_by_charge(mol_id)


def color_protein_na_for_active_molecule():
  mol_id = _active_molecule_or_status()
  if mol_id is None:
    return None
  return color_protein_na(mol_id)


def color_waters_for_active_molecule():
  mol_id = _active_molecule_or_status()
  if mol_id is None:
    return None
  return color_waters(mol_id)


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
  turn_off_backup(mol_id)
  try:
    obj_number = generic_object_with_name("chain_breaks_{mol_id}".format(mol_id=mol_id))
    generic_object_clear(obj_number)
    for ch_id in chain_ids(mol_id):
      previous_terminus_xyz = {}
      for resn in range(0, chain_n_residues(ch_id, mol_id)):
        resname = resname_from_serial_number(mol_id, ch_id, resn)
        config = _chain_break_config_for_resname(resname)
        if not config:
          continue
        atom_name, warning_threshold = config
        seqnum = seqnum_from_serial_number(mol_id, ch_id, resn)
        if seqnum == -10000:
          continue
        if (resn == 0) or (is_last_polymer_residue_sn(mol_id, ch_id, resn) == 1):
          sel_string = "//{ch_id}/{seqnum}/{atom_name}".format(ch_id=ch_id, seqnum=seqnum, atom_name=atom_name)
          make_ball_and_stick(mol_id, sel_string, 0, 0.5, 1)
        if is_term_type_mc_sn(mol_id, ch_id, resn) == 1:
          previous_xyz = _atom_xyz_or_none(mol_id, ch_id, seqnum, atom_name)
          if previous_xyz is not None:
            previous_terminus_xyz[atom_name] = previous_xyz
        if is_term_type_mn_sn(mol_id, ch_id, resn) == 1:
          previous_xyz = previous_terminus_xyz.get(atom_name)
          current_xyz = _atom_xyz_or_none(mol_id, ch_id, seqnum, atom_name)
          previous_seqnum = seqnum_from_serial_number(mol_id, ch_id, resn-1)
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
  finally:
    turn_on_backup(mol_id)
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
  user_defined_click(2,merge_chains_by_click)

#Fit all segments to map
def rigid_body_fit_segments():
  if (imol_refinement_map()==-1):
    info_dialog("You must set a refinement map!")
    return None
  mol_id = _active_molecule_or_status()
  if mol_id is None:
    return None
  previous_immediate_replacement = refinement_immediate_replacement_state()
  turn_off_backup(mol_id)
  try:
    for seg in segment_list(mol_id):
      res_start=seg[2]
      res_end=seg[3]
      ch_id=seg[1]
      set_refinement_immediate_replacement(1)
      rigid_body_refine_zone(mol_id,ch_id,res_start,res_end)
      accept_regularizement()
  finally:
    set_refinement_immediate_replacement(previous_immediate_replacement)
    turn_on_backup(mol_id)
      
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
  turn_off_backup(mol_id)
  try:
    set_refinement_immediate_replacement(1)
    rigid_body_refine_zone(
      mol_id,
      segment_context["chain_id"],
      segment_context["segment_start"],
      segment_context["segment_end"],
    )
    accept_regularizement()
  finally:
    set_refinement_immediate_replacement(previous_immediate_replacement)
    turn_on_backup(mol_id)

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
  first_in_seg = first_residue_in_seg(mol_id, ch_id, resn)
  last_in_seg = last_residue_in_seg(mol_id, ch_id, resn)
  active_atom["first_in_seg"] = first_in_seg
  active_atom["last_in_seg"] = last_in_seg
  active_atom["grow_from_n_term"] = abs(first_in_seg - resn) <= abs(last_in_seg - resn)
  return active_atom


def _run_terminal_growth_action(builder):
  growth_context = _active_terminal_growth_context_or_status()
  if not growth_context:
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
      "auto",
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
  turn_off_backup(mol_id)
  try:
    ncs_control_change_ncs_master_to_chain_id(mol_id, ch_id)
    copy_from_ncs_master_to_others(mol_id, ch_id)
  finally:
    turn_on_backup(mol_id)

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
  if len(sn_list)==1:
    sn_start=sn_list[0]
    interesting_entry=_sequence_match_list_entry(mol_id,ch_id,seq,sn_start,pattern_length)
    if not interesting_entry:
      info_dialog("Found the sequence match, but could not locate a CA/P atom to center on.")
      return None
    return _activate_interesting_entry(interesting_entry)
  elif len(sn_list)==0:
    info_dialog("Sequence not found!")
  elif len(sn_list)>1:
    for sn in sn_list:
      list_entry=_sequence_match_list_entry(mol_id,ch_id,seq,sn,pattern_length)
      if not list_entry:
        continue
      interesting_list.append(list_entry)
    if not interesting_list:
      info_dialog("Found matching sequence positions, but could not locate jump atoms for them.")
      return None
    navigable_interesting_things_gui("Matches to entered sequence",interesting_list)

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
  user_defined_click(1, grow_post_click)


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
    renumber_residue_range(mol_id,ch_id,first_res,last_res,int(offset))
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
    renumber_residue_range(mol_id,ch_id,first_res,last_res,int(offset))
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
    renumber_residue_range(mol_id,ch_id,first_res,last_res,int(offset))
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
    renumber_residue_range(mol_id,ch_id,first_res,resn,offset)
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
    renumber_residue_range(mol_id,ch_id,resn,last_res,offset)
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
  ("Common ions", [
    ("Alkalis / earths", [
      ("Li", "LI"),
      ("Na", "NA"),
      ("K", "K"),
      ("Rb", "RB"),
      ("Cs", "CS"),
      ("NH4", "NH4"),
      ("Mg", "MG"),
      ("Ca", "CA"),
      ("Sr", "SR"),
      ("Ba", "BA"),
    ]),
    ("Halides / pseudohalides", [
      ("F", "F"),
      ("Cl", "CL"),
      ("Br", "BR"),
      ("I", "IOD"),
      ("SCN", "SCN"),
      ("N3", "AZI"),
    ]),
    ("Transition metals", [
      ("Mn", "MN"),
      ("Fe", "FE"),
      ("Co", "CO"),
      ("Ni", "NI"),
      ("Cu", "CU"),
      ("Zn", "ZN"),
      ("Cd", "CD"),
      ("Hg", "HG"),
    ]),
    ("Heavier / special metals", [
      ("Pb", "PB"),
      ("Pt", "PT"),
      ("Au", "AU"),
      ("Ir", "IR"),
      ("Os", "OS"),
      ("Gd", "GD3"),
      ("La", "LA"),
      ("Sm", "SM"),
      ("Eu", "EU3"),
      ("Tb", "TB"),
      ("Ho", "HO3"),
      ("Lu", "LU"),
      ("Yb", "YB"),
      ("U", "U1"),
    ]),
  ]),
  ("Buffers / small ions", [
    ("Phosphate", "PO4"),
    ("Sulfate", "SO4"),
    ("Carbonate", "CO3"),
    ("Nitrate", "NO3"),
    ("Formate", "FMT"),
    ("Acetate", "ACT"),
    ("Citrate", "FLC"),
    ("Tartrate", "TAR"),
    ("Thiocyanate", "SCN"),
    ("Azide", "AZI"),
    ("Urea", "URE"),
  ]),
  ("Cryoprotectants / solvents", [
    ("Glycerol", "GOL"),
    ("Ethylene glycol", "EDO"),
    ("MPD", "MPD"),
    ("PEG fragment", "PEG"),
    ("DMSO", "DMS"),
    ("Ethanol", "EOH"),
    ("Isopropanol", "IPA"),
    ("Acetonitrile", "ACN"),
  ]),
  ("Common ligands", [
    ("Acetyl-CoA", "ACO"),
    ("ATP", "ATP"),
    ("ADP", "ADP"),
    ("AMP", "AMP"),
    ("GTP", "GTP"),
    ("GDP", "GDP"),
    ("NAD", "NAD"),
    ("NADP", "NAP"),
    ("FAD", "FAD"),
    ("FMN", "FMN"),
    ("HEM", "HEM"),
    ("PLP", "PLP"),
    ("CoA", "COA"),
  ]),
]


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


@lru_cache(maxsize=1)
def _ccp4_monomer_library_roots():
  roots = []
  for candidate in (
    os.environ.get("CLIBD_MON"),
    os.path.join(os.environ.get("CCP4", ""), "lib", "data", "monomers") if os.environ.get("CCP4") else "",
    "/Applications/CCP4/ccp4-9/lib/data/monomers",
  ):
    if candidate and os.path.isdir(candidate) and candidate not in roots:
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
  get_monomer_no_H(monomer_code)
  return molecule_number_list()[-1]


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

  coot.user_defined_click_py(2, make_link_from_clicks)
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

  base_chain_id, base_resno, _base_ins_code = base_spec
  new_chain_id, new_resno, _new_ins_code = new_spec
  if base_chain_id != new_chain_id:
    return 0

  status = regularize_zone(imol, base_chain_id, min(base_resno, new_resno), max(base_resno, new_resno), "")
  if status == 1:
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
      sn_max=chain_n_residues(ch_id,mol_id)
      for sn in  range(0,sn_max+1):
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
    get_monomer_no_H("ALA")
    if not model_molecule_list():
      info_dialog(f"Failed to create the starting monomer for {button_label}.")
      return 0
    mol_id=model_molecule_list()[-1]
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
    turn_off_backup(mol_id)
    try:
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
    finally:
      turn_on_backup(mol_id)
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
  chain_context = _active_chain_context_or_status()
  if not chain_context:
    return None
  current_num=chain_context["resno"]
  mol_id=chain_context["mol_id"]
  ch_id=chain_context["chain_id"]
  def renum_seg(new_num):
    new_num = _parse_int_entry(new_num, "Renumbering requires an integer residue number.")
    if new_num is None:
      return 0
    seg_count=0
    segments=segment_list(mol_id)
    last_res=last_polymer_residue(mol_id,ch_id)
    first_res=first_residue(mol_id,ch_id)
    new_seg_list=[]
    for seg in segments:
      if ch_id==seg[1]:
        new_seg_list.append(seg) 
    for seg in new_seg_list:
      seg_count=seg_count+1
      if (current_num>=seg[2]) and (current_num<=seg[3]):
        res_start=seg[2]
        res_end=seg[3]
        ch_id=seg[1]
        if res_end<last_res:
          seg_next=new_seg_list[seg_count]
        if res_start>first_res:
          seg_prev=new_seg_list[seg_count-2]
        offset=new_num-current_num
        if ((((res_start==first_res) or (res_start+offset)>seg_prev[3])) and (((res_end==last_res) or (res_end+offset)<seg_next[2]))):
          renumber_residue_range(mol_id,ch_id,res_start,res_end,int(offset))
        else:
          info_dialog("No can do, this would result in overlapping sequence numbering!")
        _refresh_extra_restraints_display(mol_id)
        return 1
    return 0
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
        imol = atom_spec_1[1]
        add_extra_bond_restraint(imol, atom_spec_1[2], atom_spec_1[3], atom_spec_1[4], atom_spec_1[5], atom_spec_1[6], atom_spec_2[2], atom_spec_2[3], atom_spec_2[4], atom_spec_2[5], atom_spec_2[6], bl, 0.035)
      user_defined_click(2, make_restr_dist)
      if continue_qm:
        user_defined_add_arbitrary_length_bond_restraint(bl)
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


def _build_custom_display_menu(submenu_display):
  """Populate the Display submenu in normal Python rather than exec'd text."""
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


def _build_custom_fit_menu(submenu_fit):
  add_simple_coot_menu_menuitem(submenu_fit, "Fit all chains to map", lambda func: rigid_fit_all_chains())
  add_simple_coot_menu_menuitem(submenu_fit, "Fit current chain to map", lambda func: rigid_fit_active_chain())
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


def build_custom_menu():
  """Create the top-level Custom menu with direct Python submenu builders."""
  coot_gui_module = _coot_gui_module()
  if coot_gui_module is None:
    return None
  if "Gio" not in globals() and hasattr(coot_gui_module, "Gio"):
    globals()["Gio"] = coot_gui_module.Gio

  menu = attach_module_menu_button("Custom")
  submenu_display = Gio.Menu.new()
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

  menu.append_submenu("Display", submenu_display)
  menu.append_submenu("Fit", submenu_fit)
  menu.append_submenu("Renumber", submenu_renumber)
  menu.append_submenu("Settings", submenu_settings)
  menu.append_submenu("Build", submenu_build)
  menu.append_submenu("Mutate", submenu_mutate)
  menu.append_submenu("Modify", submenu_modify)
  menu.append_submenu("Maps", submenu_maps)

  _build_custom_display_menu(submenu_display)
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
  return menu

if _coot_gui_module() is not None:
  try:
    build_custom_menu()
  except Exception:
    print("coot_trimmings direct menu build failed")
    traceback.print_exc()
