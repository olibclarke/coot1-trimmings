# Oli's Coot customizations.
# Copy this file into Coot's startup script directory and restart Coot.

import math
import os
import struct
import shutil
import subprocess
import sys
import traceback

import gi
import coot

if not hasattr(gi, "require_version"):
  def _require_version(_namespace, _version):
    return None
  gi.require_version = _require_version

from coot import *
from coot_utils import *
REGISTERED_KEYBINDING_CALLBACKS = []
MAP_SURFACE_DISPLAY_STATE = {}
MAP_SURFACE_OPACITY_STATE = {}
MAP_GLOBAL_VIEW_SETTINGS = {}
MAP_LOCAL_APPEARANCE_STATE = {}
SMART_COPY_TEMPLATE_IMOL = None
SMART_COPY_SOURCE_CENTRE = None
SMART_COPY_RESIDUE_NAME = None
MODEL_AMBIENT_LIGHTING_ENABLED = False
EM_REFINED_MAP_COLOUR = (0.10, 0.57, 0.95)
EM_REFINED_MAP_CONTOUR_SIGMA = 2.3
EM_TARGET_PIXEL_SIZE = 0.5

# User-adjustable startup settings.
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
STARTUP_ROTATION_CENTRE_CROSSHAIRS_SCALE = 0.15
STARTUP_ROTATION_CENTRE_CROSSHAIRS_COLOUR = (1.0, 1.0, 1.0, 1.0)
DEFAULT_NEW_HELIX_CHAIN_ID = "A"
MODEL_NORMAL_AMBIENT = (0.35, 0.35, 0.35, 1.0)
MODEL_NORMAL_DIFFUSE = (0.75, 0.75, 0.75, 1.0)
MODEL_NORMAL_SPECULAR = (0.20, 64.0)
MODEL_AMBIENT_ONLY_AMBIENT = (1.0, 1.0, 1.0, 1.0)
MODEL_AMBIENT_ONLY_DIFFUSE = (0.0, 0.0, 0.0, 1.0)
MODEL_AMBIENT_ONLY_SPECULAR = (0.0, 64.0)

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


def _optional_import(module_name):
  try:
    return __import__(module_name)
  except Exception:
    return None


coot_gui = _optional_import("coot_gui")
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
  set_user_defined_rotation_centre_crosshairs_size_scale_factor(
    STARTUP_ROTATION_CENTRE_CROSSHAIRS_SCALE
  )
  set_rotation_centre_cross_hairs_colour(
    *STARTUP_ROTATION_CENTRE_CROSSHAIRS_COLOUR
  )


def _active_residue_or_status():
  residue = active_residue()
  if residue:
    return residue
  add_status_bar_text("No active residue")
  return None


def _scrollable_map_or_status():
  map_id = scroll_wheel_map()
  if map_id != -1 and map_id in map_molecule_list():
    return map_id
  add_status_bar_text("No active map")
  return None


def _residue_is_polymer(mol_id, chain_id, resno, ins_code):
  if does_residue_exist_p(mol_id, chain_id, resno, ins_code) == 0:
    return False
  return residue_name(mol_id, chain_id, resno, ins_code) in POLYMER_RESIDUE_NAMES


def _rotation_centre_xyz():
  return [rotation_centre_position(axis) for axis in (0, 1, 2)]


def _distance_sq(point_1, point_2):
  dx = point_1[0] - point_2[0]
  dy = point_1[1] - point_2[1]
  dz = point_1[2] - point_2[2]
  return dx*dx + dy*dy + dz*dz


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


def _click_spec_imol(click_spec):
  if not isinstance(click_spec, list):
    residue = active_residue()
    return residue[0] if residue else -1
  if len(click_spec) >= 7:
    if isinstance(click_spec[0], int) and click_spec[0] in molecule_number_list():
      return click_spec[0]
    if isinstance(click_spec[1], int) and click_spec[1] in molecule_number_list():
      return click_spec[1]
  if len(click_spec) == 6 and isinstance(click_spec[0], int):
    return click_spec[0]
  residue = active_residue()
  if residue and residue[0] in molecule_number_list():
    return residue[0]
  return -1


def _click_spec_chain_id(click_spec):
  if not isinstance(click_spec, list):
    return False
  if len(click_spec) >= 7:
    return click_spec[2]
  if len(click_spec) == 6:
    return click_spec[1]
  return False


def _click_spec_res_no(click_spec):
  if not isinstance(click_spec, list):
    return False
  if len(click_spec) >= 7:
    return click_spec[3]
  if len(click_spec) == 6:
    return click_spec[2]
  return False


def _click_spec_ins_code(click_spec):
  if not isinstance(click_spec, list):
    return ""
  if len(click_spec) >= 7:
    return click_spec[4]
  if len(click_spec) == 6:
    return click_spec[3]
  return ""


def _click_spec_alt_conf(click_spec):
  if not isinstance(click_spec, list):
    return ""
  if len(click_spec) >= 7:
    return click_spec[6]
  if len(click_spec) == 6:
    return click_spec[5]
  return ""


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


if GUI_PYTHON_AVAILABLE:
  coot_toolbar_button = coot_gui.coot_toolbar_button
  attach_module_menu_button = coot_gui.attach_module_menu_button
  generic_button_dialog = coot_gui.generic_button_dialog
  generic_multiple_entries_with_check_button = (
    coot_gui.generic_multiple_entries_with_check_button
  )
  generic_single_entry = coot_gui.generic_single_entry
  add_simple_action_to_menu = coot_gui.add_simple_action_to_menu
  add_module_cryo_em_gui = coot_gui.add_module_cryo_em_gui
  add_module_refine = coot_gui.add_module_refine

  _menu_action_counter = 0

  def add_simple_coot_menu_menuitem(menu, menu_item_label, activate_function):
    global _menu_action_counter
    action_name = "coot_trimmings_action_{0}".format(_menu_action_counter)
    _menu_action_counter += 1

    def on_activate(_action, _parameter):
      return activate_function(None)

    add_simple_action_to_menu(menu, menu_item_label, action_name, on_activate)
else:
  def add_key_binding(name, key, thunk):
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

  def add_simple_coot_menu_menuitem(_menu, _menu_item_label, _activate_function):
    return None

  def add_module_cryo_em_gui():
    return None

  def add_module_refine():
    return None

if not CARBOHYDRATE_GUI_AVAILABLE:
  def add_module_carbohydrate_gui():
    return None

_apply_startup_settings()

# Gtk-dependent modules, menus, toolbar buttons, and dialog bindings are
# archived near the end of this file for clean re-enabling once the embedded
# Python Gtk layer is available again.

#place helix with prosmart alpha helix restraints and depict in rainbow
add_key_binding("Place helix here","h",
lambda: place_helix_with_restraints())

#Measure distance
add_key_binding("Measure distance","m",
lambda: start_measure_distance())

#Go to nearest map peak around the rotation centre
add_key_binding("Go to nearest density peak","b",
lambda: go_to_nearest_density_peak())

#Resample active EM map to 0.5 A/pixel and restyle
add_key_binding("Resample active EM map to 0.5 A/pixel","C",
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

add_key_binding("Map plus 0.5 sigma","|",
lambda: step_current_map_coarse_up())

add_key_binding("Map minus 0.5 sigma","_",
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
# def toggle_all_models_ambient_lighting():
#   global MODEL_AMBIENT_LIGHTING_ENABLED
#   if MODEL_AMBIENT_LIGHTING_ENABLED:
#     ambient = MODEL_NORMAL_AMBIENT
#     diffuse = MODEL_NORMAL_DIFFUSE
#     specular = MODEL_NORMAL_SPECULAR
#     status_message = "Set all models to normal lighting"
#   else:
#     ambient = MODEL_AMBIENT_ONLY_AMBIENT
#     diffuse = MODEL_AMBIENT_ONLY_DIFFUSE
#     specular = MODEL_AMBIENT_ONLY_SPECULAR
#     status_message = "Set all models to ambient lighting"
#   for imol in model_molecule_number_list():
#     set_model_material_ambient(imol, ambient[0], ambient[1], ambient[2], ambient[3])
#     set_model_material_diffuse(imol, diffuse[0], diffuse[1], diffuse[2], diffuse[3])
#     set_model_material_specular(imol, specular[0], specular[1])
#   MODEL_AMBIENT_LIGHTING_ENABLED = not MODEL_AMBIENT_LIGHTING_ENABLED
#   add_status_bar_text(status_message)


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


def smart_paste_copied_non_polymer_residue():
  if SMART_COPY_TEMPLATE_IMOL not in molecule_number_list():
    add_status_bar_text("No copied ligand, ion, or water available")
    return None
  residue = active_residue()
  if residue:
    target_mol_id = residue[0]
  else:
    target_mol_id = go_to_atom_molecule_number()
    if target_mol_id not in model_molecule_list():
      add_status_bar_text("No active model for smart paste")
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
  try:
    translate_molecule_by(paste_imol, dx, dy, dz)
    merge_result = merge_molecules_py([paste_imol], target_mol_id)
    if not merge_result or merge_result[0] != 1:
      add_status_bar_text("Unable to merge pasted residue into active model")
      return None
  finally:
    if paste_imol in molecule_number_list():
      close_molecule(paste_imol)
  add_status_bar_text(
    "Pasted {resname} at pointer".format(
      resname=(SMART_COPY_RESIDUE_NAME or "copied residue")
    )
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

  residue_serial_entries = all_residues_with_serial_numbers_py(mol_id) or []
  sequence_index_by_residue = {}
  next_sequence_index_by_chain = {}
  ordered_residue_specs_by_chain = {}
  for residue_entry in residue_serial_entries:
    if not residue_entry or len(residue_entry) < 4:
      continue
    chain_id, resno, ins_code = residue_entry[1:4]
    sequence_index_by_residue[(chain_id, resno, ins_code)] = next_sequence_index_by_chain.get(chain_id, 0)
    next_sequence_index_by_chain[chain_id] = next_sequence_index_by_chain.get(chain_id, 0) + 1
    ordered_residue_specs_by_chain.setdefault(chain_id, []).append((chain_id, resno, ins_code))

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
      if not nearby_spec or len(nearby_spec) < 3:
        continue
      residue_spec_2 = tuple(nearby_spec[:3])
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
  if imol_refinement_map()==-1:
    add_status_bar_text("You need to set a refinement map")
    return None
  residue = _active_residue_or_status()
  if not residue:
    return None
  ch_id = residue[1]
  resno = residue[2]
  altloc = residue[5]
  previous_immediate_replacement = refinement_immediate_replacement_state()
  try:
    set_refinement_immediate_replacement(1)
    refine_zone(mol_id, ch_id, resno, resno, altloc)
    accept_regularizement()
  finally:
    set_refinement_immediate_replacement(previous_immediate_replacement)


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


def increase_map_radius():
  current_radius = get_map_radius()
  set_map_radius(current_radius + 2.0)
  set_map_radius_em(current_radius + 2.0)


def decrease_map_radius():
  current_radius = get_map_radius()
  new_radius = max(2.0, current_radius - 2.0)
  set_map_radius(new_radius)
  set_map_radius_em(new_radius)


def active_map_surface_displayed():
  map_id = _scrollable_map_or_status()
  if map_id is None:
    return False
  return (
    map_is_displayed(map_id)
    and MAP_SURFACE_DISPLAY_STATE.get(map_id, False)
  )


def increase_active_map_surface_opacity():
  if not active_map_surface_displayed():
    add_status_bar_text("Active map is not in solid-surface mode")
    return None
  map_id = scroll_wheel_map()
  current_opacity = MAP_SURFACE_OPACITY_STATE.get(
    map_id, get_solid_density_surface_opacity(map_id)
  )
  new_opacity = min(1.0, current_opacity + 0.05)
  MAP_SURFACE_OPACITY_STATE[map_id] = new_opacity
  set_solid_density_surface_opacity(map_id, new_opacity)


def decrease_active_map_surface_opacity():
  if not active_map_surface_displayed():
    add_status_bar_text("Active map is not in solid-surface mode")
    return None
  map_id = scroll_wheel_map()
  current_opacity = MAP_SURFACE_OPACITY_STATE.get(
    map_id, get_solid_density_surface_opacity(map_id)
  )
  new_opacity = max(0.0, current_opacity - 0.05)
  MAP_SURFACE_OPACITY_STATE[map_id] = new_opacity
  set_solid_density_surface_opacity(map_id, new_opacity)


def display_only_active_map():
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
  try:
    mol_id_active=active_residue()[0]
  except:
    mol_id_active=model_molecule_list()[0]
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
    
def step_map_coarse_up(mol_id):
  current_level=get_contour_level_in_sigma(mol_id)
  if (current_level >= 0.5) and (current_level <= 10.0):
    new_level=current_level+0.5
  elif (current_level<0.5):
    new_level=0.5
  elif (current_level>10.0):
    new_level=10.0
  set_contour_level_in_sigma(mol_id,new_level)

def step_map_coarse_down(mol_id):
  current_level=get_contour_level_in_sigma(mol_id)
  if (current_level >= 0.5) and (current_level <= 10.0):
    new_level=current_level-0.5
  elif (current_level<0.5):
    new_level=0.5
  elif (current_level>10.0):
    new_level=10.0
  set_contour_level_in_sigma(mol_id,new_level)


def set_current_map_sigma(level):
  map_id = _scrollable_map_or_status()
  if map_id is None:
    return None
  set_contour_level_in_sigma(map_id, level)


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


def start_measure_distance():
  add_status_bar_text("Pick 2 atoms to measure distance")
  do_distance_define()


def _map_looks_em_like(map_id):
  if not is_valid_map_molecule(map_id):
    return False
  if map_is_difference_map(map_id):
    return False
  try:
    cell = map_cell(map_id)
  except Exception:
    return False
  if not cell or len(cell) < 6:
    return False
  alpha, beta, gamma = cell[3], cell[4], cell[5]
  return (
    abs(alpha - 90.0) < 0.5
    and abs(beta - 90.0) < 0.5
    and abs(gamma - 90.0) < 0.5
  )


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
  if map_is_difference_map(map_id):
    add_status_bar_text(status_message)
    return None
  style_resampled_em_map(map_id)
  set_contour_level_in_sigma(map_id, contour_level_sigma)
  set_map_displayed(map_id, 1)
  set_scroll_wheel_map(map_id)
  set_scrollable_map(map_id)
  add_status_bar_text(status_message)
  return map_id


def resample_active_map_for_em_half_angstrom(force=False):
  map_id = _scrollable_map_or_status()
  if map_id is None:
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
  resample_factor, current_pixel_size = _resample_factor_for_target_pixel_size(
    map_id, EM_TARGET_PIXEL_SIZE
  )
  if resample_factor is None:
    return _restyle_active_map_without_resampling(
      map_id,
      contour_level_sigma,
      "Grid spacing could not be determined",
    )
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
    set_map_radius(fallback_radius)
    return
  set_map_radius(saved["map_radius"])
  set_clipping_front(saved["clipping_front"])
  set_clipping_back(saved["clipping_back"])


def _set_global_view_extent(map_id):
  current_radius = max(get_map_radius(), 1.0)
  global_radius = max(_map_global_extent_radius(map_id), current_radius)
  front = get_clipping_plane_front()
  back = get_clipping_plane_back()
  scale = global_radius / current_radius
  set_map_radius(global_radius)
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
  residue = _active_residue_or_status()
  if not residue:
    return None
  mol_id=residue[0]
  ch_id=residue[1]
  resn=residue[2]
  atom_name=residue[4]
  sn=get_sn_from_resno(mol_id,ch_id,resn)
  next_sn=sn+1
  next_res=seqnum_from_serial_number(mol_id,"%s"%(ch_id),next_sn)
  set_go_to_atom_molecule(mol_id)
  if (next_res!=-10000 and is_protein_chain_p(mol_id,ch_id)==1):
    set_go_to_atom_chain_residue_atom_name(ch_id,next_res,atom_name)
  elif (next_res!=-10000 and is_nucleotide_chain_p(mol_id,ch_id)==1):
    set_go_to_atom_chain_residue_atom_name(ch_id,next_res,atom_name)

#Go to previous residue in current polymer chain.
def prev_res():
  residue = _active_residue_or_status()
  if not residue:
    return None
  mol_id=residue[0]
  ch_id=residue[1]
  resn=residue[2]
  atom_name=residue[4]
  sn=get_sn_from_resno(mol_id,ch_id,resn)
  if (sn>=1):
    sn=sn-1
  prev_res=seqnum_from_serial_number(mol_id,"%s"%(ch_id),sn)
  set_go_to_atom_molecule(mol_id)
  if (prev_res!=-10000 and is_protein_chain_p(mol_id,ch_id)==1):
    set_go_to_atom_chain_residue_atom_name(ch_id,prev_res,atom_name)
  elif (prev_res!=-10000 and is_nucleotide_chain_p(mol_id,ch_id)==1):
    set_go_to_atom_chain_residue_atom_name(ch_id,prev_res,atom_name)
  
def sequence_context():
  mol_id=active_residue()[0]
  ch_id=active_residue()[1]
  resnum=active_residue()[2]
  ins_code=active_residue()[3]
  resname=residue_name(mol_id,ch_id,resnum,ins_code)
  def get_aa_code(resnum):
    mol_id=active_residue()[0]
    ch_id=active_residue()[1]
    ins_code=active_residue()[3]
    if residue_name(mol_id,ch_id,resnum,ins_code):
      aa_code=three_letter_code2single_letter(residue_name(mol_id,ch_id,resnum,ins_code))
      if (len(residue_name(mol_id,ch_id,resnum,ins_code))==1):
        aa_code=residue_name(mol_id,ch_id,resnum,ins_code)
      if (residue_name(mol_id,ch_id,resnum,ins_code)!="ALA") and (aa_code=="A") and (len(residue_name(mol_id,ch_id,resnum,ins_code))!=1):
        aa_code="X"
    else:
      aa_code="-"
    return aa_code
  current_res=get_aa_code(resnum)
  minus_context=""
  for i in range(1,10):
    aa_code=str(get_aa_code(resnum-i))
    minus_context=aa_code + minus_context
  plus_context=""
  for i in range(1,10):
    aa_code=str(get_aa_code(resnum+i))
    plus_context=plus_context + aa_code
  final_string="Residue:  " + resname + "  " + str(resnum) +"  Sequence context: ..." + minus_context + "[" + current_res + "]" + plus_context + "..."
  info_dialog(final_string)
  
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

#Open displayed models and maps in Chimera, set map levels and view
def open_in_chimera():
  if find_exe("chimera"): #If chimera exists do stuff, else raise an error
    import subprocess
    pwd=os.getcwd() #Dir from which coot was launched
    cofr=rotation_centre() #3-membered list [x,y,z]
    make_directory_maybe("coot-chimera") #Put coot droppings here
    coot_chimera_path=pwd+"/coot-chimera/"
    check_path=coot_chimera_path+"chimera_launcher.cmd" #Chimera run script that will be written later
    check_path2=coot_chimera_path+"matrix.txt" #Orientation matrix for chimera input
    view_number=add_view_here("tmp")
#     initial_zoom=zoom_factor() #Coot's zoom factor. 100 is the initial state, zooming out gets larger, in gets smaller.
#     reset_view()
    zoom=zoom_factor()
#     go_to_view_number(view_number,1)
#     scale_factor=initial_zoom/base_zoom
#     chimera_scale=base_zoom/initial_zoom
    chimera_scale=100.0/zoom #Approx - breaks down when initial views upon loading mol in Coot or Chimera are different
    coot_slab_thickness=(0.071*zoom)+2.8612 #Empirically estimated
    chimera_clip_val=coot_slab_thickness/2.0 #Dist from cofr to clip plane (symmetric)
    map_radius=get_map_radius() #We'll use this to export a map_radius sized fragment
    if os.path.isfile(check_path): #if chimera script exists, delete it.
      os.remove(check_path)
    if os.path.isfile(check_path2): #same for matrix.txt
      os.remove(check_path2)
    map_list=map_molecule_list()
    mol_list=model_molecule_list()
    disp_mol_list=[]
    disp_map_list=[]
    for mol_id in model_molecule_list():
      if mol_is_displayed(mol_id):
        disp_mol_list.append(mol_id) 
    for map_id in map_molecule_list():
      if map_is_displayed(mol_id):
        disp_map_list.append(map_id)   
    print(("disp model list",disp_mol_list))
    print(("disp map list",disp_map_list))
    set_graphics_window_size(1000,1000)
    matrix_list=view_matrix() #9-membered list of rotation matrix elements.
    with open(check_path2,"a") as matrix_file: #Write orientation matrix (specify first pdb as model, then use matrixcopy to apply to all others). Translation vector in last column is 0,0,0 because we will take care of that separately.
      matrix_file.write("Model {model_id}.0\n".format(model_id=disp_mol_list[0]))
      matrix_file.write("\t {a} {b} {c} 0.0\n".format(a=matrix_list[0],b=matrix_list[1],c=matrix_list[2]))
      matrix_file.write("\t {a} {b} {c} 0.0\n".format(a=matrix_list[3],b=matrix_list[4],c=matrix_list[5]))
      matrix_file.write("\t {a} {b} {c} 0.0\n".format(a=matrix_list[6],b=matrix_list[7],c=matrix_list[8]))
    with open(check_path,"a") as cmd_file: #Start writing stuff to chimera launch script.
      for mol_id in model_molecule_list():
        if mol_is_displayed(mol_id):
          file_name=coot_chimera_path+"mol_{mol_id}.pdb".format(mol_id=mol_id)
          write_pdb_file(mol_id,file_name)
          path_to_file=file_name
          model_id=mol_id
          cmd_file.write("open #{model_id} {mol}\n".format(model_id=model_id,mol=path_to_file)) 
          cmd_file.write("color gold #{model_id}; color byhet #{model_id};  sel #{model_id}; namesel tmp; ~ribbon tmp; ~disp tmp; sel tmp&@CA&protein; repr stick sel; disp sel; sel tmp&~protein; repr stick sel; disp sel;  setattr m stickScale 1.0 tmp;  sel side chain/base.without CA/C1'&tmp; repr stick sel; disp sel; setattr b radius 0.1 sel; color byhet tmp; sel @CA|@CB; namesel tmp2; sel tmp&tmp2; repr stick sel; setattr b radius 0.1 sel; sel @CD,N&:pro&tmp; disp sel; repr stick sel; setattr b radius 0.1 sel; ~sel; disp ~protein\n".format(model_id=model_id)) 
          cmd_file.write("matrixset matrix.txt\n")
          cmd_file.write("matrixcopy #{mol_id0} #{mol_id}\n".format(mol_id0=disp_mol_list[0],mol_id=model_id))         
      map_id0=mol_id
      for map_id in map_molecule_list():
        if map_is_displayed(map_id):
          map_level=get_contour_level_absolute(map_id)
          if map_is_difference_map(map_id):
            model_id=map_id0+map_id #make sure that the assigned model number of each map is unique and does not overlap with those of pdbs.
            if map_id==0:
              model_id=model_id+1
            file_name=coot_chimera_path+"diff_map_{model_id}.mrc".format(model_id=model_id)
            path_to_file=file_name
            export_map_fragment(map_id,cofr[0],cofr[1],cofr[2],map_radius,file_name)
            cmd_file.write("open #{model_id} {diff_map} \n".format(model_id=model_id,diff_map=path_to_file))
            cmd_file.write("volume #{model_id} capfaces false style mesh meshlighting false squaremesh false color \"#08882eefa222\" step 1 ;  sop cap off ;  set depthCue ;  set dcStart 0.2 ;  set dcEnd 1 ;  background solid white ; set showcofr ;  cofr view ;  clip on; volume #{model_id} level -{map_level} color #da1200000000 level {map_level} color #0000bda00000 \n".format(model_id=model_id,map_level=map_level))
            cmd_file.write("matrixset matrix.txt\n")
            cmd_file.write("matrixcopy #{mol_id0} #{model_id}\n".format(mol_id0=disp_mol_list[0],model_id=model_id))
          else:
            model_id=map_id0+map_id
            if map_id==0:
              model_id=model_id+1
            file_name=coot_chimera_path+"map_{model_id}.mrc".format(model_id=model_id)
            export_map_fragment(map_id,cofr[0],cofr[1],cofr[2],map_radius,file_name)
            path_to_file=file_name
            cmd_file.write("open #{model_id} {map} \n".format(model_id=model_id,map=path_to_file))
            cmd_file.write("volume #{model_id} capfaces false style mesh meshlighting false squaremesh false color \"#08882eefa222\" step 1 ;  sop cap off ;  set depthCue ;  set dcStart 0.2 ;  set dcEnd 1 ;  background solid white ; set showcofr ;  cofr view ;  clip on; volume #{model_id} level {map_level} \n".format(model_id=model_id,map_level=map_level))
            cmd_file.write("matrixset matrix.txt\n")
            cmd_file.write("matrixcopy #{mol_id0} #{model_id}\n".format(mol_id0=disp_mol_list[0],model_id=model_id))
      cmd_file.write("~sel; set projection orthographic; clip off; cofr {x},{y},{z} coordinatesystem #{mol_id0}; ac mc; center sel; cofr view; cofr fixed; clip hither {clipval} fromCenter true; clip yon -{clipval} fromCenter true; cofr view; clip on; ~set showcofr; del sel; scale {chimera_scale}; windowsize 1000 1000; scene ca_and_sidechains save; disp; scene all_atoms save; scene ca_and_sidechains reset; windowsize 1000 1000".format(x=cofr[0],y=cofr[1],z=cofr[2],mol_id0=disp_mol_list[0],chimera_scale=chimera_scale,clipval=chimera_clip_val))
      #line above applies translation and scale obtained from coot rotation center and zoom_factor
    chimera_exe=find_exe("chimera")
    info_dialog("Opening in Chimera...\nOrientation should be right but you will probably need to \nadjust the scale and clipping to get the same view as in Coot. \n If you can't see anything, try zooming out.")
    subprocess.Popen(["chimera",check_path])
  else: 
    info_dialog("Sorry, you need UCSF Chimera installed and accessible from the terminal for this to work!")
#This works okay, though adjustment of zoom and clipping still not ideal.
#would it be possible to add a subprocess.communicate() to use Coot as a controller for Chimera?
#Hmmm
#Doesn't work when attempting to open second set of maps loaded, when first set is undisplayed... Fixed! (I think)


#Colour active segment
def colour_active_segment():
  mol_id=active_residue()[0]
  segments=segment_list(mol_id)
  res_here=active_residue()[2]
  ins_code=active_residue()[3]
  ch_id=active_residue()[1]
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
        res_color_spec=[([ch_id,res,""],segment_colour)]
        colour_list=colour_list+res_color_spec
    else:
      res_start=seg[2]
      res_end=seg[3]
      ch_id_here=seg[1]
      for res in range(res_start,res_end+1):
        blank_color_spec=[([ch_id_here,res,""],blank_colour)]
        blank_list=blank_list+blank_color_spec
  clear_user_defined_atom_colours(mol_id)
  set_user_defined_atom_colour_by_residue_py(mol_id,colour_list)
  set_user_defined_atom_colour_by_residue_py(mol_id,blank_list)
  graphics_to_user_defined_atom_colours_representation(mol_id)


def color_emringer_outliers(mol_id,map_id):
  if find_exe("phenix.emringer"):
    import subprocess
    import sys
    mmtbx_path=find_exe("phenix")[:-16]+"modules/cctbx_project/"
    sys.path.insert(0,mmtbx_path)
    import mmtbx
    import pickle as pickle
    pwd=os.getcwd()
    model_name=molecule_name(mol_id)
    map_name=molecule_name(map_id)
    make_directory_maybe("coot-emringer")
    coot_emringer_path=pwd+"/coot-emringer/"
    output_file_name=coot_emringer_path+"mol_{mol_id}_cablam_output.txt".format(mol_id=mol_id)
    p=subprocess.Popen("phenix.emringer {model_name} {map_name}".format(model_name=model_name,map_name=map_name),shell=True)
    p.communicate()
    emringer_outlier_list=[]
    emringer_outlier_color=30
    emringer_outlier_pkl=pwd+"/"+molecule_name_stub(mol_id,2)+"_emringer_plots/Outliers.pkl"
    outlier_string=str(pickle.load(open(emringer_outlier_pkl,"rb")))
    with open(output_file_name,"a") as outlier_file: 
      outlier_file.write(outlier_string)
    with open(output_file_name) as f:
      emringer_output = [x.strip('\n') for x in f.readlines()] #make a list of lines in cablam output file, stripping newlines
      for line in emringer_output:
        line_list=line.split()
        if len(line_list)>=5:
          ch_id=str(line_list[2]) # string.split() defaults to splitting by spaces and making into a list
          print(("ch_id:",ch_id))
          resid=int(line_list[1]) #If second column has trailing letters (as it will if there is an insertion code) then strip them
          print(("resid",resid))
          ins_id=""
          emringer_outlier_color_spec=[([ch_id,resid,ins_id],emringer_outlier_color)]
          emringer_outlier_list=emringer_outlier_list+emringer_outlier_color_spec
      print(("outlier_list",emringer_outlier_list))
      try:
        set_user_defined_atom_colour_by_residue_py(mol_id,emringer_outlier_list)
        graphics_to_user_defined_atom_colours_representation(mol_id)
      except NameError:
        info_dialog("You need a newer Coot - custom coloring is only in r6174 and later, sorry.")
        pass
  else:
    info_dialog("Sorry, you need phenix.cablam_validate, sed and awk installed and accessible from the terminal for this to work!")

def color_by_cablam2(mol_id):
  if find_exe("phenix.cablam_validate") and find_exe("awk") and find_exe("sed"):
    import subprocess
    pwd=os.getcwd() #dir where coot was launched
    make_directory_maybe("coot-cablam") #Put coot droppings here
    coot_cablam_path=pwd+"/coot-cablam/"
    file_name=molecule_name(mol_id)
    file_name_output=coot_cablam_path+"mol_{mol_id}_cablam_output.txt".format(mol_id=mol_id) #this file will have info on cablam outliers
#     write_pdb_file(mol_id,file_name) # write a copy of the active mol to the cablam dir
    p=subprocess.Popen("phenix.cablam_validate {file_name} outlier_cutoff=0.01 | tail -n +2 | awk '{{FS=\":\"}} {{print $1,$2,$5}}' | sed 's/CaBLAM Outlier/1/g' | sed 's/CaBLAM Disfavored/2/g' | sed 's/CA Geom Outlier/3/g' | sed 's/try alpha helix/4/g' | sed 's/try beta sheet/5/g' | sed 's/try three-ten/6/g' | sed 's/[0123456789-]/ &/' > {file_name_output}".format(file_name=file_name,file_name_output=file_name_output),shell=True)
    p.communicate() #wait for cablam to finish
    cablam_outlier_list=[]
    cablam_beta_list=[]
    cablam_alpha_list=[]
    cablam_3_10_list=[]
    cablam_disfavored_list=[]
    cablam_outlier_colour=30
    cablam_disfavored_colour=27
    cablam_alpha_colour=10
    cablam_beta_colour=39
    cablam_3_10_colour=15
    cablam_ca_geom_colour=10
    blank_colour=0
    outlier_flag=0
    with open(file_name_output) as f:
      cablam_output = [x.strip('\n') for x in f.readlines()] #make a list of lines in cablam output file, stripping newlines
      for line in cablam_output:
        line_list=line.split()
        if len(line_list)>=3:
          ch_id=line_list[0] # string.split() defaults to splitting by spaces and making into a list
          resid=int(line_list[1].rstrip('abcdefghjklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ')) #If second column has trailing letters (as it will if there is an insertion code) then strip them
          ins_id=str(line_list[1].lstrip('0123456789')) #If second column has trailing characters, these correspond to the insertion code. Strip leading numbers.
          resname=line_list[2]
          if len(line_list)>3:
            outlier_flag=int(line_list[3]) #1 is Cablam outlier 2 is cablam disfavored, and 3 is ca geom outlier
            if outlier_flag==1: #CaBLAM outliers
              cablam_outlier_colour_spec=[([ch_id,resid,ins_id],cablam_outlier_colour)]
              cablam_outlier_list=cablam_outlier_list+cablam_outlier_colour_spec
            elif outlier_flag==2: #CaBLAM disfavored
              cablam_outlier_colour_spec=[([ch_id,resid,ins_id],cablam_disfavored_colour)]
              cablam_outlier_list=cablam_outlier_list+cablam_outlier_colour_spec
            elif outlier_flag==3: #CaBLAM disfavored
              cablam_outlier_colour_spec=[([ch_id,resid,ins_id],cablam_ca_geom_colour)]
              cablam_outlier_list=cablam_outlier_list+cablam_outlier_colour_spec
            elif outlier_flag==4: #alpha
              cablam_outlier_colour_spec=[([ch_id,resid,ins_id],cablam_alpha_colour)]
              cablam_outlier_list=cablam_outlier_list+cablam_outlier_colour_spec
            elif outlier_flag==5: #beta
              cablam_outlier_colour_spec=[([ch_id,resid,ins_id],cablam_beta_colour)]
              cablam_outlier_list=cablam_outlier_list+cablam_outlier_colour_spec
            elif outlier_flag==6: #3-10 helix
              cablam_outlier_colour_spec=[([ch_id,resid,ins_id],cablam_3_10_colour)]
              cablam_outlier_list=cablam_outlier_list+cablam_outlier_colour_spec
            else:
              cablam_outlier_colour_spec=[([ch_id,resid,ins_id],blank_colour)]
              cablam_outlier_list=cablam_outlier_list+cablam_outlier_colour_spec
          else:
            cablam_outlier_colour_spec=[([ch_id,resid,ins_id],blank_colour)]
            cablam_outlier_list=cablam_outlier_list+cablam_outlier_colour_spec 
#     os.remove(file_name)
    os.remove(file_name_output)
    try:
      set_user_defined_atom_colour_by_residue_py(mol_id,cablam_outlier_list)
      graphics_to_user_defined_atom_colours_representation(mol_id)
      info_dialog("CaBLAM coloring scheme (predicted SS and outliers): \n  \n Teal = alpha \n \n Green = 3-10 \n \n Purple = beta \n \n Yellow = Ca geometry outlier \n \n Orange = CaBLAM disfavored (5% cutoff) \n \n Red = CaBLAM outlier (1% cutoff)")
    except NameError:
      info_dialog("You need a newer Coot - custom coloring is only in r6174 and later, sorry.")
      pass
  else:
    info_dialog("Sorry, you need phenix.cablam_validate, sed and awk installed and accessible from the terminal for this to work!")

def color_by_rama(mol_id):
  if find_exe("phenix.ramalyze") and find_exe("awk") and find_exe("sed"):
    import subprocess
    pwd=os.getcwd() #dir where coot was launched
    make_directory_maybe("coot-ramalyze") #Put coot droppings here
    coot_ramalyze_path=pwd+"/coot-ramalyze/"
    file_name=molecule_name(mol_id)
    file_name_output=coot_ramalyze_path+"mol_{mol_id}_ramalyze_output.txt".format(mol_id=mol_id) #this file will have info on ramalyze outliers
#     write_pdb_file(mol_id,file_name) # write a copy of the active mol to the ramalyze dir
    p=subprocess.Popen("phenix.ramalyze {file_name} outliers_only=true rama_potential=emsley | awk '{{FS=\":\"}} {{print $1}}' | sed 's/[0123456789-]/ &/' | awk '{{if (NF==3) print}}' > {file_name_output}".format(file_name=file_name,file_name_output=file_name_output),shell=True)
    p.communicate() #wait for ramalyze to finish
    ramalyze_outlier_list=[]
    ramalyze_outlier_colour=34
    blank_res_list=[]
    blank_colour=0
    for ch_id in chain_ids(mol_id):
      sn_max=chain_n_residues(ch_id,mol_id)
      for sn in range(0,sn_max+1):
        resn=seqnum_from_serial_number(mol_id,ch_id,sn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
        residue_to_color=[([ch_id,resn,ins_id],blank_colour)]
        blank_res_list=blank_res_list+residue_to_color
    with open(file_name_output) as f:
      ramalyze_output = [x.strip('\n') for x in f.readlines()] #make a list of lines in ramalyze output file, stripping newlines
      for line in ramalyze_output:
        line_list=line.split()
        ch_id=line_list[0] # string.split() defaults to splitting by spaces and making into a list
        resid=int(line_list[1].rstrip('abcdefghjklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ')) #If second column has trailing letters (as it will if there is an insertion code) then strip them
        ins_id=str(line_list[1].lstrip('0123456789')) #If second column has trailing characters, these correspond to the insertion code. Strip leading numbers.
        resname=line_list[2]
        ramalyze_outlier_colour_spec=[([ch_id,resid,ins_id],ramalyze_outlier_colour)]
        ramalyze_outlier_list=ramalyze_outlier_list+ramalyze_outlier_colour_spec
#     os.remove(file_name)
    os.remove(file_name_output)
    try:
      set_user_defined_atom_colour_by_residue_py(mol_id,blank_res_list)
      set_user_defined_atom_colour_by_residue_py(mol_id,ramalyze_outlier_list)
      graphics_to_user_defined_atom_colours_representation(mol_id)
    except NameError:
      info_dialog("You need a newer Coot - custom coloring is only in r6174 and later, sorry.")
      pass
  else:
    info_dialog("Sorry, you need phenix.ramalyze, sed and awk installed and accessible from the terminal for this to work!")

def color_by_cc(mol_id):
  if find_exe("phenix.model_map_cc") and find_exe("awk") and find_exe("sed") and find_exe("grep") and imol_refinement_map()!=-1 and space_group(imol_refinement_map())=="P 1":
    import subprocess
    pwd=os.getcwd() #dir where coot was launched
    make_directory_maybe("coot-cc") #Put coot droppings here
    coot_cc_path=pwd+"/coot-cc/"
    pdb_name=molecule_name(mol_id)
    map_id=int(imol_refinement_map())
    map_name=molecule_name(map_id)
#     write_pdb_file(mol_id,pdb_name) # write a copy of the active mol to the cc dir
    if not (map_name.endswith(".ccp4") or map_name.endswith(".mrc") or map_name.endswith(".map")):
      map_name=coot_cc_path+"mol_{mol_id}.ccp4".format(mol_id=mol_id)
      export_map(map_id,map_name)
    file_name_output=coot_cc_path+"mol_{mol_id}_cc_output.txt".format(mol_id=mol_id) #this file will have info on cc outliers
    p=subprocess.Popen("phenix.model_map_cc {pdb_name} {map_name} resolution=4.0 | awk '{{if (NF==7) print}}' | grep chain > {file_name_output}".format(pdb_name=pdb_name,map_name=map_name,file_name_output=file_name_output),shell=True)
    p.communicate() #wait for cc to finish
    cc_list=[]
    with open(file_name_output) as f:
      cc_output = [x.strip('\n') for x in f.readlines()] #make a list of lines in cc output file, stripping newlines
      for line in cc_output:
        line_list=line.split()
        ch_id=line_list[2] # string.split() defaults to splitting by spaces and making into a list
        resid=int(line_list[4].rstrip('abcdefghjklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ')) #If second column has trailing letters (as it will if there is an insertion code) then strip them
        ins_id=str(line_list[4].lstrip('0123456789')) #If second column has trailing characters, these correspond to the insertion code. Strip leading numbers.
        cc=float(line_list[6])
        if cc<0.0:
          cc=0.0
        cc_colour=int((1.0-cc)*31+2)
        print(("cc=",cc,"color=",cc_colour))
        cc_colour_spec=[([ch_id,resid,ins_id],cc_colour)]
        cc_list=cc_list+cc_colour_spec
    try:
      set_user_defined_atom_colour_by_residue_py(mol_id,cc_list)
      graphics_to_user_defined_atom_colours_representation(mol_id)
    except NameError:
      info_dialog("You need a newer Coot - custom coloring is only in r6188 and later, sorry.")
      pass
  else:
    info_dialog("Sorry, you need phenix.model_map_cc, sed and awk installed and accessible from the terminal for this to work! Also, needs to be P1 map right now, sorry.")
    

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
cycle_rep_flag={0:0}
def cycle_rep_up(mol_id,flag):
  global cycle_rep_flag
  cycle_rep_flag[mol_id]=flag
  if cycle_rep_flag[mol_id]==0:
    graphics_to_ca_plus_ligands_representation(mol_id)
    cycle_rep_flag[mol_id]=1
  elif cycle_rep_flag[mol_id]==1:
    graphics_to_ca_plus_ligands_and_sidechains_representation(mol_id)
    cycle_rep_flag[mol_id]=2
  elif cycle_rep_flag[mol_id]==2:
    graphics_to_ca_plus_ligands_sec_struct_representation(mol_id)
    cycle_rep_flag[mol_id]=3
  elif cycle_rep_flag[mol_id]==3:
    graphics_to_rainbow_representation(mol_id)
    cycle_rep_flag[mol_id]=4
  elif cycle_rep_flag[mol_id]==4:
    graphics_to_bonds_representation(mol_id)
    cycle_rep_flag[mol_id]=5
  elif cycle_rep_flag[mol_id]==5:
    graphics_to_b_factor_representation(mol_id)
    cycle_rep_flag[mol_id]=6
  elif cycle_rep_flag[mol_id]==6:
    try:
      graphics_to_user_defined_atom_colours_representation(mol_id)
      cycle_rep_flag[mol_id]=7
    except NameError:
      cycle_rep_flag[mol_id]=0
  elif cycle_rep_flag[mol_id]==7:
    try:
      graphics_to_user_defined_atom_colours_all_atoms_representation(mol_id)
      cycle_rep_flag[mol_id]=0
    except NameError:
      cycle_rep_flag[mol_id]=0

# def add_partial_water():
#   place_typed_atom_at_pointer("Water")
#   refine_active_residue()
#   mol_id=active_residue()[0]
#   ch_id=active_residue()[1]
#   resno=active_residue()[2]
#   ins_code=active_residue()[3]
#   atom_name=active_residue()[4]
#   alt_conf=""
#   set_atom_attribute(mol_id,ch_id,resno,ins_code,atom_name,alt_conf,"occ",0.5)
# add_key_binding("Add partial water","w",
# lambda: add_partial_water())

# def delete_alt_confs():
#   mol_id=active_residue()[0]
#   for ch_id in chain_ids(mol_id):
#     for resn in (first_residue(mol_id,ch_id),last_residue(mol_id,ch_id)+1):
#      residue_inf=residue_info(mol_id,ch_id,resn,"")      

    
def cycle_rep_down(mol_id,flag):
  global cycle_rep_flag
  cycle_rep_flag[mol_id]=flag
  if cycle_rep_flag[mol_id]==4:
    graphics_to_ca_plus_ligands_sec_struct_representation(mol_id)
    cycle_rep_flag[mol_id]=3
  elif cycle_rep_flag[mol_id]==3:
    graphics_to_ca_plus_ligands_and_sidechains_representation(mol_id)
    cycle_rep_flag[mol_id]=2
  elif cycle_rep_flag[mol_id]==2:
    graphics_to_ca_plus_ligands_representation(mol_id)
    cycle_rep_flag[mol_id]=1
  elif cycle_rep_flag[mol_id]==1:
    try:
      graphics_to_user_defined_atom_colours_all_atoms_representation(mol_id)
      cycle_rep_flag[mol_id]=0
    except NameError:
      cycle_rep_flag[mol_id]=7
  elif cycle_rep_flag[mol_id]==0:
    try:
      graphics_to_user_defined_atom_colours_representation(mol_id)
      cycle_rep_flag[mol_id]=7
    except NameError:
      cycle_rep_flag[mol_id]=7
  elif cycle_rep_flag[mol_id]==7:
    graphics_to_b_factor_representation(mol_id)
    cycle_rep_flag[mol_id]=6
  elif cycle_rep_flag[mol_id]==6:
    graphics_to_bonds_representation(mol_id)
    cycle_rep_flag[mol_id]=5
  elif cycle_rep_flag[mol_id]==5:
    graphics_to_rainbow_representation(mol_id)
    cycle_rep_flag[mol_id]=4


def cycle_rep_up_current():
  residue = _active_residue_or_status()
  if not residue:
    return None
  mol_id = residue[0]
  cycle_rep_up(mol_id, cycle_rep_flag.get(mol_id, 0))


def cycle_rep_down_current():
  residue = _active_residue_or_status()
  if not residue:
    return None
  mol_id = residue[0]
  cycle_rep_down(mol_id, cycle_rep_flag.get(mol_id, 0))


#Refine triple (Paul)
def key_binding_refine_triple():
    active_atom = active_residue()
    if not active_atom:
       print("No active atom")
    else:
       imol = active_atom[0]
       residue_spec = atom_spec_to_residue_spec(active_atom)
       N_terminal_residue_spec = [residue_spec_to_chain_id(residue_spec),
                                  residue_spec_to_res_no(residue_spec)-1,
                                  residue_spec_to_ins_code(residue_spec)]
       C_terminal_residue_spec = [residue_spec_to_chain_id(residue_spec),
                                  residue_spec_to_res_no(residue_spec)+1,
                                  residue_spec_to_ins_code(residue_spec)]
       spec_list = [N_terminal_residue_spec, residue_spec, C_terminal_residue_spec]
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
  update_go_to_atom_from_current_position()
  if rotamer_number>=n_rots:
    rotamer_number=0
    set_residue_to_rotamer_number(mol_id,ch_id,res_here,ins_code,alt_conf,rotamer_number)
  else:
    rotamer_number=rotamer_number+1
    set_residue_to_rotamer_number(mol_id,ch_id,res_here,ins_code,alt_conf,rotamer_number)
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
  sn=0
  resn_to_match=active_residue()[2]
  mol_id=active_residue()[0]
  ch_id=active_residue()[1]
  current_resn=seqnum_from_serial_number(mol_id,"%s"%(ch_id),sn)
  while (current_resn!=resn_to_match):
    sn=sn+1
    current_resn=seqnum_from_serial_number(mol_id,"%s"%(ch_id),sn)
  return sn
    
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
    active_spec = [chain_id, residue_no, ins_code]
    if is_polymer_residue_spec(molecule_id, tuple(active_spec)):
      return [chain_id, residue_no, ins_code]

    for search_radius in (4.0, 6.0, 8.0, 10.0):
      nearby_specs = residues_near_residue(molecule_id, active_spec, search_radius)
      polymer_specs = [spec for spec in nearby_specs if spec and len(spec) == 3 and is_polymer_residue_spec(molecule_id, tuple(spec))]
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
      if nearby_residue and len(nearby_residue) == 3:
        nearby_spec = tuple(nearby_residue)
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
#Deletes active chain
def delete_chain():
  active_chain_id=active_residue()[1]
  active_mol_id=active_residue()[0]
  turn_off_backup(active_mol_id)
  while (is_polymer(active_mol_id,active_chain_id)==1) or (
  is_solvent_chain_p(active_mol_id,active_chain_id)!=-1):
    first_res=first_residue(active_mol_id,active_chain_id)
    last_res=last_residue(active_mol_id,active_chain_id)
    delete_residue_range(active_mol_id,active_chain_id,first_res,last_res)
  turn_on_backup(active_mol_id)

#Fits all polymer chains to map
def rigid_fit_all_chains():
  mol_id=active_residue()[0]
  turn_off_backup(mol_id)
  for ch_id in chain_ids(mol_id): #Rigid body refine each chain
    if is_polymer(mol_id,ch_id)==1:
      rigid_body_refine_by_atom_selection(mol_id, "//%s//"%(ch_id))
      accept_regularizement()
  turn_on_backup(mol_id)
  
#Fits active chain to map
def rigid_fit_active_chain():
  mol_id=active_residue()[0]
  ch_id=active_residue()[1]
  rigid_body_refine_by_atom_selection(mol_id, "//%s//"%(ch_id))
  
#Copies active chain
def copy_active_chain():
  mol_id=active_residue()[0]
  ch_id=active_residue()[1]
  new_molecule_by_atom_selection(mol_id, "//%s//"%(ch_id))
  
#Cuts active chain
def cut_active_chain():
  mol_id=active_residue()[0]
  ch_id=active_residue()[1]
  new_molecule_by_atom_selection(mol_id, "//%s//"%(ch_id))
  turn_off_backup(mol_id)
  while (is_polymer(mol_id,ch_id)==1) or (is_solvent_chain_p(mol_id,ch_id)!=-1):
    first_res=first_residue(mol_id,ch_id)
    last_res=last_residue(mol_id,ch_id)
    delete_residue_range(mol_id,ch_id,first_res,last_res)
  turn_on_backup(mol_id)

#Faster rigid body fit (not relevant for smaller ranges, but much faster for larger ranges)
def fast_rigid_fit(res_start,res_end,ch_id,mol_id):
  if (mol_id in model_molecule_number_list()) and (ch_id in chain_ids(mol_id)):
    turn_off_backup(mol_id)
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
    turn_on_backup(mol_id)
  else:
    info_dialog("The specified chain or molecule does not exist!")

#Copy active segment
def copy_active_segment():
  mol_id=active_residue()[0]
  segments=segment_list(mol_id)
  res_here=active_residue()[2]
  ch_id=active_residue()[1]
  for seg in segments:
    if (res_here>=seg[2]) and (res_here<=seg[3]) and (ch_id==seg[1]):
      res_start=seg[2]
      res_end=seg[3]
      ch_id=seg[1]
      new_molecule_by_atom_selection(mol_id, "//{ch_id}/{res_start}-{res_end}/".format(ch_id=ch_id,res_start=res_start,res_end=res_end))

#Cut active segment
def cut_active_segment():
  mol_id=active_residue()[0]
  segments=segment_list(mol_id)
  res_here=active_residue()[2]
  ch_id=active_residue()[1]
  for seg in segments:
    if (res_here>=seg[2]) and (res_here<=seg[3]) and (ch_id==seg[1]):
      res_start=seg[2]
      res_end=seg[3]
      ch_id=seg[1]
      new_molecule_by_atom_selection(mol_id, "//{ch_id}/{res_start}-{res_end}/".format(ch_id=ch_id,res_start=res_start,res_end=res_end))
      delete_residue_range(mol_id,ch_id,res_start,res_end)

#Delete active segment
def delete_active_segment():
  mol_id=active_residue()[0]
  segments=segment_list(mol_id)
  res_here=active_residue()[2]
  ch_id=active_residue()[1]
  for seg in segments:
    if (res_here>=seg[2]) and (res_here<=seg[3]) and (ch_id==seg[1]):
      res_start=seg[2]
      res_end=seg[3]
      ch_id=seg[1]
      delete_residue_range(mol_id,ch_id,res_start,res_end)


#Jiggle-fits active chain to map
def jiggle_fit_active_chain():
  if (imol_refinement_map()==-1):
    info_dialog("You must set a refinement map!")
  else:
    mol_id=active_residue()[0]
    ch_id=active_residue()[1]
    fit_chain_to_map_by_random_jiggle(mol_id,ch_id,1000,0.1)

#Jiggle-fit active chain to B-smoothed map
def jiggle_fit_active_chain_smooth():
  if (imol_refinement_map()==-1):
    info_dialog("You must set a refinement map!")
  else:
    sharpen(imol_refinement_map(),200)
    mol_id=active_residue()[0]
    ch_id=active_residue()[1]
    fit_chain_to_map_by_random_jiggle(mol_id,ch_id,1000,0.1)
    sharpen(imol_refinement_map(),0)
    
#Jiggle-fits active chain to map (more thorough)
def jiggle_fit_active_chain_slow():
  if (imol_refinement_map()==-1):
    info_dialog("You must set a refinement map!")
  else:
    mol_id=active_residue()[0]
    ch_id=active_residue()[1]
    score_string=""
    for scale_fac in range(1,1000,100):
      scale_fac2=scale_fac/10.0
      score=fit_chain_to_map_by_random_jiggle(mol_id,ch_id,100,scale_fac2)
      score_string=score_string+"scale:"+str(scale_fac2)+"score:"+str(score)
    print(score_string)
    
#Fits all polymer chains to map
def jiggle_fit_all_chains():
  if (imol_refinement_map()==-1):
    info_dialog("You must set a refinement map!")
  else:
    mol_id=active_residue()[0]
    turn_off_backup(mol_id)
    for ch_id in chain_ids(mol_id): #Rigid body refine each chain
      if is_polymer(mol_id,ch_id)==1:
        fit_chain_to_map_by_random_jiggle(mol_id,ch_id,1000,0.1)
        accept_regularizement()
    turn_on_backup(mol_id)

#Jiggle-fit current molecule to map
def jiggle_fit_active_mol():
  if (imol_refinement_map()==-1):
    info_dialog("You must set a refinement map!")
  else:
    mol_id=active_residue()[0]
    fit_molecule_to_map_by_random_jiggle(mol_id,1000,0.1)
    
#Clear labels and distances
def clear_distances_and_labels():
  remove_all_atom_labels()
  clear_measure_distances()
  
#Delete hydrogens from active molecule
def delete_h_active():
  mol_id=active_residue()[0]
  delete_hydrogens(mol_id)
  
#click the start and end point, then fit the gap between them with polyala
def fit_polyala_gui():
  def fit_polyala(res1,res2):
    length=abs(res1[3]-res2[3])-1
    loop_seq=length*"A"
    fit_gap(res1[1],res1[2],res1[3],res2[3],loop_seq,1)
  user_defined_click(2,fit_polyala)
  
# Try to rebuild with db_mainchain after fit_gap?
# def fit_polyala_gui2():
#   def fit_polyala(res1,res2):
#     length=abs(res1[3]-res2[3])-1
#     loop_seq=length*"A"
#     fit_gap(res1[1],res1[2],res1[3],res2[3],loop_seq,1)
#     db_mainchain?
#   user_defined_click(2,fit_polyala)

#Rebuild backbone in selected zone
def rebuild_backbone_wrapper():
  def rebuild_backbone(res1,res2):
    if res1[1]==res2[1] and res1[2]==res2[2]: #if residues in same mol and chain
      mol_id=res1[1]
      ch_id=res1[2]
      resid1=res1[3]
      resid2=res2[3]
      if resid1!=resid2:
        if resid2<resid1:
          resid1, resid2 = resid2, resid1
        new_mol_id=db_mainchain(mol_id,ch_id,resid1,resid2,"forwards")
        accept_regularizement()
        res1_rsr=first_residue(new_mol_id,ch_id)
        res2_rsr=last_residue(new_mol_id,ch_id)
        #need to mutate each residue to appropriate sidechain and kill sidechain
        #get target seq with aa_code=three_letter_code2single_letter(residue_name(mol_id,ch_id,resnum,ins_code))
        #mutate_residue_range
        #delete_sidechain_range
        target_seq=""
        for res in range(res1_rsr,res2_rsr+1):
          aa_code=three_letter_code2single_letter(residue_name(new_mol_id,ch_id,res,""))
          target_seq=target_seq+aa_code
        print(("target seq:",target_seq))
        mutate_residue_range(new_mol_id,ch_id,res1_rsr,res2_rsr,target_seq)
        delete_sidechain_range(new_mol_id,ch_id,res1_rsr,res2_rsr)
        for res in range(res1_rsr,res2_rsr+1):
          if residue_name(new_mol_id,ch_id,res,"")=="PRO":
            target_seq="P"
            mutate_residue_range(new_mol_id,ch_id,res,res,target_seq)
        #cut out orginal region and merge in new fragment?
        refine_zone(new_mol_id,ch_id,res1_rsr,res2_rsr,"")
        accept_regularizement()
      else:
        info_dialog("Sorry, you need at least 2 residues in a zone!")
    else:
      info_dialog("Sorry, residues must be in same mol and chain!")
  user_defined_click(2,rebuild_backbone)

def rebuild_backbone_reverse_wrapper():
  def rebuild_backbone_reverse(res1):
    mol_id=res1[1]
    ch_id=res1[2]
    resid1=res1[3]
    segments=segment_list(mol_id)
    for seg in segments:
      if (resid1>=seg[2]) and (resid1<=seg[3]) and (ch_id==seg[1]):
        res_start=seg[2]
        res_end=seg[3]
        ch_id=seg[1]
        turn_off_backup(mol_id)
        if res_start!=res_end:
          reverse_direction_of_fragment(mol_id,ch_id,res_start)
          new_mol_id=db_mainchain(mol_id,ch_id,res_start,res_end,"forwards")
        accept_regularizement()
        res1_rsr=first_residue(new_mol_id,ch_id)
        res2_rsr=last_residue(new_mol_id,ch_id)
        #need to mutate each residue to appropriate sidechain and kill sidechain
        #get target seq with aa_code=three_letter_code2single_letter(residue_name(mol_id,ch_id,resnum,ins_code))
        #mutate_residue_range
        #delete_sidechain_range
        target_seq=""
        for res in range(res1_rsr,res2_rsr+1):
          aa_code=three_letter_code2single_letter(residue_name(new_mol_id,ch_id,res,""))
          target_seq=target_seq+aa_code
        print(("target seq:",target_seq))
        mutate_residue_range(new_mol_id,ch_id,res1_rsr,res2_rsr,target_seq)
        delete_sidechain_range(new_mol_id,ch_id,res1_rsr,res2_rsr)
        for res in range(res1_rsr,res2_rsr+1):
          if residue_name(new_mol_id,ch_id,res,"")=="PRO":
            target_seq="P"
            mutate_residue_range(new_mol_id,ch_id,res,res,target_seq)
        #cut out orginal region and merge in new fragment?
        delete_residue_range(mol_id,ch_id,res_start,res_end)
        merge_molecules([new_mol_id],mol_id)
        turn_off_backup(mol_id)
  user_defined_click(1,rebuild_backbone_reverse)

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

  
  
#refine range, plus everything in contact with it
def refine_residues_sphere_click():
  def refine_residues_sphere_click_a(res1,res2):
    mol_id_1=_click_spec_imol(res1)
    mol_id_2=_click_spec_imol(res2)
    ch_id_1=_click_spec_chain_id(res1)
    ch_id_2=_click_spec_chain_id(res2)
    resno_1=_click_spec_res_no(res1)
    resno_2=_click_spec_res_no(res2)
    ins_code_1=_click_spec_ins_code(res1)
    ins_code_2=_click_spec_ins_code(res2)
    if (mol_id_1==mol_id_2) and (ch_id_1==ch_id_2) and (imol_refinement_map()!=-1):
      if resno_1>=resno_2:
        resno_1,resno_2=resno_2,resno_1 #Swap res numbers if wrong way round
      res_list=[]
      nearby_residues=[]
      for resn in range(resno_1,resno_2+1):
        res_triple=[ch_id_1,resn,ins_code_1] #Does not account for varying ins code. Should be able to fix using residues_matching_criteria(), e.g.:
        #residues_matching_criteria(0, lambda chain_id,resno,ins_code,serial: True), except substitute a test func for true
        #that evaluates to true if resno, mol id and ch_id match, and returns ins_code (third item in output triple)
        res_list.append(res_triple) #append rather than adding, bc making list of lists
      for resn in range(resno_1+1,resno_2-1):
        res_triple=[ch_id_1,resn,ins_code_1] #Does not account for varying ins code. Should be able to fix using residues_matching_criteria(), e.g.:
        new_nearby_residues=residues_near_residue(mol_id_1,res_triple,5)
        nearby_residues.extend(new_nearby_residues)
      if len(res_list)<=2:
        for resn in range(resno_1,resno_2+1):
          res_triple=[ch_id_1,resn,ins_code_1] #Does not account for varying ins code. Should be able to fix using residues_matching_criteria(), e.g.:
          new_nearby_residues=residues_near_residue(mol_id_1,res_triple,5)
          nearby_residues.extend(new_nearby_residues)
      res_list=res_list+nearby_residues
      print(res_list)
      res_list=[list(x) for x in set(tuple(x) for x in res_list)] #Uniquifies list of lists
      print(("hopefully unique:",res_list))
      refine_residues(mol_id_1,res_list)
      print(res_list)
    else:
      info_dialog("Sorry, start and end residues must be in same mol and chain! And you need to have a refinement map.")
  user_defined_click(2,refine_residues_sphere_click_a)

#Python version of sphere refine from here:
#https://www.mail-archive.com/coot@jiscmail.ac.uk/msg01463.html
def sphere_refine_active(radius):
    active_atom = active_residue()
    if (not active_atom):
        add_status_bar_text("No active residue")
    else:
        imol      = active_atom[0]
        chain_id  = active_atom[1]
        res_no    = active_atom[2]
        ins_code  = active_atom[3]
        atom_name = active_atom[4]
        alt_conf  = active_atom[5]
        centred_residue = active_atom[1:4]
        other_residues = residues_near_residue(imol, centred_residue, radius)
        all_residues = [centred_residue]
        if isinstance(other_residues, list):
            all_residues += other_residues
        print("imol: %s residues: %s" %(imol, all_residues))
        refine_residues(imol, all_residues)

def stepped_sphere_refine(mol_id,ch_id):
  turn_off_backup(mol_id)
  set_refinement_immediate_replacement(1)
  valid_resnames=['A','C','T','G','U','ALA','UNK','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','ILE','LEU','LYS','MET','MSE','PHE','PRO','SER','THR','TRP','TYR','VAL']
  if is_polymer(mol_id,ch_id): 
    first_res=first_residue(mol_id,ch_id)
    last_res=last_residue(mol_id,ch_id)
    if is_protein_chain_p(mol_id,ch_id):
      for res in range(first_res,last_res+1):
        res_exist_flag=does_residue_exist_p(mol_id,ch_id,res,"")
        resname=residue_name(mol_id,ch_id,res,"")
        if (res_exist_flag==1) and (resname in valid_resnames):
          set_go_to_atom_chain_residue_atom_name(ch_id,res," CA ")
          sphere_refine_active(7)
          accept_regularizement()
    elif is_nucleotide_chain_p(mol_id,ch_id):
      for res in range(first_res,last_res+1):
        res_exist_flag=does_residue_exist_p(mol_id,ch_id,res,"")
        resname=residue_name(mol_id,ch_id,res,"")
        if (res_exist_flag==1) and (resname in valid_resnames):
          set_go_to_atom_chain_residue_atom_name(ch_id,res," P ")
          sphere_refine_active(7)
          accept_regularizement()
  set_refinement_immediate_replacement(0)
  turn_on_backup(mol_id)
  info_dialog("Refinement finished - all done!")

  
#Stepped v of cylinder refine:
#get mol from active residue
#for ch_id in mol:
#cylinder_refine four res segments, with two res overlap
#ignore breaks?
#what about waters?


#Copy fragment (click start and end)
def copy_frag_by_click():
  def copy_frag(res1,res2):
    mol_id_1=res1[1]
    mol_id_2=res2[1]
    ch_id_1=res1[2]
    ch_id_2=res2[2]
    resno_1=res1[3]
    resno_2=res2[3]
    if (mol_id_1==mol_id_2) and (ch_id_1==ch_id_2):
      if resno_1>resno_2:
        atom_sel="//%s/%s-%s" %(ch_id_1,resno_2,resno_1)
        new_molecule_by_atom_selection(mol_id_1,atom_sel)
      elif resno_2>resno_1:
        atom_sel="//%s/%s-%s" %(ch_id_1,resno_1,resno_2)
        new_molecule_by_atom_selection(mol_id_1,atom_sel)
      else:
        atom_sel="//%s/%s" %(ch_id_1,resno_2)
        new_molecule_by_atom_selection(mol_id_1,atom_sel)
    else:
      info_dialog("Start and end residues must be in the same molecule and chain!")
  user_defined_click(2,copy_frag)
  
#Cut fragment (click start and end)
def cut_frag_by_click():
  def cut_frag(res1,res2):
    mol_id_1=res1[1]
    mol_id_2=res2[1]
    ch_id_1=res1[2]
    ch_id_2=res2[2]
    resno_1=res1[3]
    resno_2=res2[3]
    if (mol_id_1==mol_id_2) and (ch_id_1==ch_id_2):
      turn_off_backup(mol_id_1)
      if resno_1>resno_2:
        atom_sel="//%s/%s-%s" %(ch_id_1,resno_2,resno_1)
        new_molecule_by_atom_selection(mol_id_1,atom_sel)
        delete_residue_range(mol_id_1,ch_id_1,resno_2,resno_1)
      elif resno_2>resno_1:
        atom_sel="//%s/%s-%s" %(ch_id_1,resno_1,resno_2)
        new_molecule_by_atom_selection(mol_id_1,atom_sel)
        delete_residue_range(mol_id_1,ch_id_1,resno_1,resno_2)
      else:
        atom_sel="//%s/%s" %(ch_id_1,resno_2)
        new_molecule_by_atom_selection(mol_id_1,atom_sel)
        delete_residue_range(mol_id_1,ch_id_1,resno_2)
      turn_on_backup(mol_id_1)
    else:
      info_dialog("Start and end residues must be in the same molecule and chain!")
  user_defined_click(2,cut_frag)


#Delete sidechain range (click start and end)
def delete_sidechain_range_by_click_a():
  def delete_sidechain_range_by_click_b(res1,res2):
    if (res1[1]!=res2[1]) or (res1[2]!=res2[2]) or (res1[3]==res2[3]):
      info_dialog("Start and end points must be in the same mol and chain!")
    else: 
      if (res1[3] > res2[3]):
        res_start=res2[3]-1
        res_end=res1[3]+1
      else:
        res_start=res1[3]-1
        res_end=res2[3]+1
      mol_id=res1[1]
      ch_id=res1[2]
      turn_off_backup(mol_id)
      delete_sidechain_range(mol_id,ch_id,res_start,res_end)
      turn_on_backup(mol_id)
  user_defined_click(2,delete_sidechain_range_by_click_b)
  
#Mutate range to poly-unk
def mutate_residue_range_by_click_a():
  def mutate_residue_range_by_click_b(res1,res2):
    if (res1[1]!=res2[1]) or (res1[2]!=res2[2]) or (res1[3]==res2[3]):
      info_dialog("Start and end points must be in the same mol and chain!")
    else:
      if (res1[3] > res2[3]):
        res_start=res2[3]
        res_end=res1[3]
        n=res_end-res_start+1
      else:
        res_start=res1[3]
        res_end=res2[3]
        n=res_end-res_start+1
      mol_id=res1[1]
      ch_id=res1[2]
      target_seq=n*"A"
      turn_off_backup(mol_id)
      mutate_residue_range(mol_id,ch_id,res_start,res_end,target_seq)
      for n in range(res_start,res_end+1):
        set_residue_name(mol_id,ch_id,n,"","UNK")
      turn_on_backup(mol_id)
  user_defined_click(2,mutate_residue_range_by_click_b)

#Mutate range to polyala
def mutate_residue_range_by_click_ala_a():
  def mutate_residue_range_by_click_ala_b(res1,res2):
    if (res1[1]!=res2[1]) or (res1[2]!=res2[2]) or (res1[3]==res2[3]):
      info_dialog("Start and end points must be in the same mol and chain!")
    else:
      if (res1[3] > res2[3]):
        res_start=res2[3]
        res_end=res1[3]
        n=res_end-res_start+1
      else:
        res_start=res1[3]
        res_end=res2[3]
        n=res_end-res_start+1
      mol_id=res1[1]
      ch_id=res1[2]
      target_seq=n*"A"
      turn_off_backup(mol_id)
      mutate_residue_range(mol_id,ch_id,res_start,res_end,target_seq)
      for n in range(res_start,res_end+1):
        set_residue_name(mol_id,ch_id,n,"","ALA")
      turn_on_backup(mol_id)
  user_defined_click(2,mutate_residue_range_by_click_ala_b)
  
#Merge two fragments
def merge_fragments():
  def merge_2_fragments(res1,res2):
    mol_daughter=[res2[1]]
    mol_ref=res1[1]
    merge_molecules(mol_daughter,mol_ref)
    toggle_display_mol(mol_ref)
    toggle_display_mol(mol_ref)
  user_defined_click(2,merge_2_fragments)
  
#Force addition of residue - useful when
#Coot says "No acceptable position found"
# but density is clear.
def force_add_terminal_residue():
  def force_addition(res1):
    mol_id=res1[1]
    ch_id=res1[2]
    res_no=res1[3]
    res_type="auto"
    add_terminal_residue_using_phi_psi(mol_id,ch_id,res_no,
    res_type,-57.82,-47)
    if residue_exists_qm(mol_id,ch_id,res_no+1,res1[4]):
      set_b_factor_residue_range(mol_id,ch_id,res_no+1,res_no+1,default_new_atoms_b_factor())
    elif residue_exists_qm(mol_id,ch_id,res_no-1,res1[4]):
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
       print("No active atom")
    else:
       imol = active_atom[0]
       residue_spec = atom_spec_to_residue_spec(active_atom)
       print(('spin_N {} {} {}'.format(imol, residue_spec, 120)))
       spin_N_py(imol, residue_spec, 120)

#This works, but should alter to be more flexible. ideally, have one key to tweak phi, one to tweak psi.
#Will need two global cycle variables - residue_phi_cycle and residue_psi_cycle - as well as two variables to keep the current value of phi and psi.
residue_phi_cycle=0
residue_psi_cycle=0
current_phi=-60
current_psi=-50
#should chamnge this to incorporate measurement of starting phi/psi using get_torsion (and maybe setting using set_torsion?)
#get_torsion(0,["A",2393,""," C  ",""], ["A",2394,"", " N  ", ""], ["A", 2394, "", " CA ", ""], ["A", 2394, "", " C  ",""])
#set_torsion(imol, chain_id, res_no, ins_code_alt_conf, atom_name_1,atom_name_2, atom_name_3, atom_name_4)
def cycle_residue_phi():
  global residue_phi_cycle
  global current_phi
  global current_psi
  res_type="auto"
  mol_id=active_residue()[0]
  ch_id=active_residue()[1]
  resn=active_residue()[2]
  atom_name=active_residue()[4]
  set_go_to_atom_molecule(mol_id)
  ins_code=""
  first_in_seg=first_residue_in_seg(mol_id,ch_id,resn)
  last_in_seg=last_residue_in_seg(mol_id,ch_id,resn)
  delta_first=abs(first_in_seg-resn)
  delta_last=abs(last_in_seg-resn)
  set_new_atom_b_fac_to_mean()
  if delta_first<=delta_last:
    set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
    if (residue_phi_cycle==0):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=-180
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_phi_cycle=1
    elif (residue_phi_cycle==1):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=-150
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_phi_cycle=2
    elif (residue_phi_cycle==2):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=-120
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_phi_cycle=3
    elif (residue_phi_cycle==3):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=-90
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_phi_cycle=4
    elif (residue_phi_cycle==4):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=-60
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_phi_cycle=5
    elif (residue_phi_cycle==5):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=-30
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_phi_cycle=6
    elif (residue_phi_cycle==6):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=0
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_phi_cycle=7
    elif (residue_phi_cycle==7):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=30
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_phi_cycle=8
    elif (residue_phi_cycle==8):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=60
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_phi_cycle=9
    elif (residue_phi_cycle==9):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=90
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_phi_cycle=10
    elif (residue_phi_cycle==10):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=120
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_phi_cycle=0
  else:
    set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
    if (residue_phi_cycle==0):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=-180
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_phi_cycle=1
    elif (residue_phi_cycle==1):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=-150
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_phi_cycle=2
    elif (residue_phi_cycle==2):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=-120
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_phi_cycle=3
    elif (residue_phi_cycle==3):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=-90
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_phi_cycle=4
    elif (residue_phi_cycle==4):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=-60
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_phi_cycle=5
    elif (residue_phi_cycle==5):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=-30
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_phi_cycle=6
    elif (residue_phi_cycle==6):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=0
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_phi_cycle=7
    elif (residue_phi_cycle==7):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=30
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_phi_cycle=8
    elif (residue_phi_cycle==8):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=60
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_phi_cycle=9
    elif (residue_phi_cycle==9):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=90
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_phi_cycle=10
    elif (residue_phi_cycle==10):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=120
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_phi_cycle=0
  current_phi=phi
def cycle_residue_psi():
  global residue_psi_cycle
  global current_phi
  global current_psi
  res_type="auto"
  mol_id=active_residue()[0]
  ch_id=active_residue()[1]
  resn=active_residue()[2]
  atom_name=active_residue()[4]
  set_go_to_atom_molecule(mol_id)
  ins_code=""
  first_in_seg=first_residue_in_seg(mol_id,ch_id,resn)
  last_in_seg=last_residue_in_seg(mol_id,ch_id,resn)
  delta_first=abs(first_in_seg-resn)
  delta_last=abs(last_in_seg-resn)
  set_new_atom_b_fac_to_mean()
  if delta_first<=delta_last:
    set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
    if (residue_psi_cycle==0):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=current_phi
      psi=-180
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_psi_cycle=1
    elif (residue_psi_cycle==1):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=current_phi
      psi=-150
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_psi_cycle=2
    elif (residue_psi_cycle==2):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=current_phi
      psi=-120
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_psi_cycle=3
    elif (residue_psi_cycle==3):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=current_phi
      psi=-90
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_psi_cycle=4
    elif (residue_psi_cycle==4):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=current_phi
      psi=-60
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_psi_cycle=5
    elif (residue_psi_cycle==5):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=current_phi
      psi=-30
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_psi_cycle=6
    elif (residue_psi_cycle==6):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=current_phi
      psi=0
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_psi_cycle=7
    elif (residue_psi_cycle==7):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=current_phi
      psi=30
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_psi_cycle=8
    elif (residue_psi_cycle==8):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=current_phi
      psi=60
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_psi_cycle=9
    elif (residue_psi_cycle==9):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=current_phi
      psi=90
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_psi_cycle=10
    elif (residue_psi_cycle==10):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=current_phi
      psi=120
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_psi_cycle=0
  else:
    set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
    if (residue_psi_cycle==0):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=current_phi
      psi=-180
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_psi_cycle=1
    elif (residue_psi_cycle==1):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=current_phi
      psi=-150
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_psi_cycle=2
    elif (residue_psi_cycle==2):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=current_phi
      psi=-120
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_psi_cycle=3
    elif (residue_psi_cycle==3):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=current_phi
      psi=-90
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_psi_cycle=4
    elif (residue_psi_cycle==4):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=current_phi
      psi=-60
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_psi_cycle=5
    elif (residue_psi_cycle==5):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=current_phi
      psi=-30
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_psi_cycle=6
    elif (residue_psi_cycle==6):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=current_phi
      psi=0
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_psi_cycle=7
    elif (residue_psi_cycle==7):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=current_phi
      psi=30
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_psi_cycle=8
    elif (residue_psi_cycle==8):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=current_phi
      psi=60
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_psi_cycle=9
    elif (residue_psi_cycle==9):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=current_phi
      psi=90
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_psi_cycle=10
    elif (residue_psi_cycle==10):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=current_phi
      psi=120
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_psi_cycle=0
  current_psi=psi


#Mutate mets to MSE
def mutate_all_mets_to_mse():
  mol_id=active_residue()[0]
  turn_off_backup(mol_id)
  for ch_id in chain_ids(mol_id):
    for resn in range(0,chain_n_residues(ch_id,mol_id)):
      if (resname_from_serial_number(mol_id,ch_id,resn)=="MET"):
        seqnum=seqnum_from_serial_number(mol_id,ch_id,resn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,resn))
        delete_residue_with_full_spec(mol_id,1,ch_id,seqnum,ins_id,"B")
        delete_residue_with_full_spec(mol_id,1,ch_id,seqnum,ins_id,"C")
        mutate(mol_id,ch_id,seqnum,ins_id,"MSE")
        map_id=imol_refinement_map()
        if (map_id!=-1):
          auto_fit_best_rotamer(seqnum,"",ins_id,ch_id,mol_id,map_id,1,0.01)
  turn_on_backup(mol_id)
  
#Mutate MSEs to MET
def mutate_all_mse_to_met():
  mol_id=active_residue()[0]
  turn_off_backup(mol_id)
  for ch_id in chain_ids(mol_id):
    for resn in range(0,chain_n_residues(ch_id,mol_id)):
      if (resname_from_serial_number(mol_id,ch_id,resn)=="MSE"):
        seqnum=seqnum_from_serial_number(mol_id,ch_id,resn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,resn))
        delete_residue_with_full_spec(mol_id,1,ch_id,seqnum,ins_id,"B")
        delete_residue_with_full_spec(mol_id,1,ch_id,seqnum,ins_id,"C")
        mutate(mol_id,ch_id,seqnum,ins_id,"MET")
        map_id=imol_refinement_map()
        if (map_id!=-1):
          auto_fit_best_rotamer(seqnum,"",ins_id,ch_id,mol_id,map_id,1,0.01)
  turn_on_backup(mol_id)
  
#Shorten loop by one residue
def shorten_loop():
  active_atom=active_residue()
  mol_id=active_atom[0]
  ch_id=active_atom[1]
  resn=active_atom[2]
  ins_code=active_atom[3]
  delete_residue(mol_id,ch_id,resn,ins_code)
  first_res=seqnum_from_serial_number(mol_id,"%s"%(ch_id),0)
  renumber_residue_range(mol_id,ch_id,first_res,resn,1)
  delete_all_extra_restraints(mol_id)
  set_show_extra_restraints(mol_id,0)
  set_show_extra_restraints(mol_id,1)
  r1=resn-1
  r2=resn+2
  set_refinement_immediate_replacement(1)
  refine_zone(mol_id,ch_id,r1,r2,"")
  accept_regularizement()
  set_refinement_immediate_replacement(0)

#Lengthen loop by one residue
def lengthen_loop():
  active_atom=active_residue()
  mol_id=active_atom[0]
  ch_id=active_atom[1]
  resn=active_atom[2]
  first_res=seqnum_from_serial_number(mol_id,"%s"%(ch_id),0)
  renumber_residue_range(mol_id,ch_id,first_res,resn,-1)
  delete_all_extra_restraints(mol_id)
  set_show_extra_restraints(mol_id,0)
  set_show_extra_restraints(mol_id,1)
  r1=resn-1
  r2=resn
  set_refinement_immediate_replacement(1)
  fit_gap(mol_id,ch_id,r1,r2,"A",1)
  accept_regularizement()
  set_refinement_immediate_replacement(0)

#Get fractional coordinates of active atom. Useful when inspecting heavy atom sites.
def get_fract_coords():
  a=active_residue()
  x_cart=atom_specs(a[0],a[1],a[2],a[3],a[4],a[5])[3]
  y_cart=atom_specs(a[0],a[1],a[2],a[3],a[4],a[5])[4]
  z_cart=atom_specs(a[0],a[1],a[2],a[3],a[4],a[5])[5]
  mol_id=active_residue()[0]
  cell_a=cell(mol_id)[0]
  cell_b=cell(mol_id)[1]
  cell_c=cell(mol_id)[2]
  cell_alpha=math.radians(cell(mol_id)[3])
  cell_beta=math.radians(cell(mol_id)[4])
  cell_gamma=math.radians(cell(mol_id)[5])
  cos_alpha_star=(math.cos(cell_beta)*math.cos(cell_gamma)-math.cos(cell_alpha))/(math.sin(cell_beta)*math.sin(cell_gamma))
  sin_alpha_star=math.sqrt(1-cos_alpha_star**2)
  z_fract=z_cart/(cell_c*math.sin(cell_beta)*sin_alpha_star)
  y_fract=(y_cart-(-1*cell_c*math.sin(cell_beta)*cos_alpha_star)*z_fract)/(cell_b*math.sin(cell_gamma))
  x_fract=(x_cart-(cell_b*math.cos(cell_gamma))*y_fract-(cell_c*math.cos(cell_beta))*z_fract)/cell_a
  x_y_z_string=str("("+str("%.3f" % x_fract)+","+str("%.3f" % y_fract)+","+str("%.3f" % z_fract)+")")
  info_dialog("The fractional coordinates of the active atom are: %s"%(x_y_z_string))
  
#Go to center of scroll wheel map.
def goto_center_of_map():
  if (scroll_wheel_map()==-1):
    info_dialog("You need a map!")
  else:
    mol_id=scroll_wheel_map()
    a=float(cell(mol_id)[0])
    b=float(cell(mol_id)[1])
    c=float(cell(mol_id)[2])
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
  mol_id=active_residue()[0]
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
  try:
    set_user_defined_atom_colour_by_residue_py(mol_id,rotamer_outlier_list)
    set_user_defined_atom_colour_by_residue_py(mol_id,missing_atoms_list)
    graphics_to_user_defined_atom_colours_representation(mol_id)
  except NameError:
    info_dialog("You need a newer Coot - custom coloring is only in r6174 and later, sorry.")
    pass


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
  try:
    set_user_defined_atom_colour_by_residue_py(mol_id,polar_res_list)
    set_user_defined_atom_colour_by_residue_py(mol_id,hphob_res_list)
    set_user_defined_atom_colour_by_residue_py(mol_id,gly_res_list)
    set_user_defined_atom_colour_by_residue_py(mol_id,pro_res_list)
    set_user_defined_atom_colour_by_residue_py(mol_id,blank_res_list)
    graphics_to_user_defined_atom_colours_representation(mol_id)
  except NameError:
    info_dialog("You need a newer Coot - custom coloring is only in r6174 and later, sorry.")
    pass
  
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
  try:
    set_user_defined_atom_colour_by_residue_py(mol_id,pos_res_list)
    set_user_defined_atom_colour_by_residue_py(mol_id,neg_res_list)
    set_user_defined_atom_colour_by_residue_py(mol_id,blank_res_list)
    graphics_to_user_defined_atom_colours_representation(mol_id)
  except NameError:
    info_dialog("You need a newer Coot - custom coloring is only in r6174 and later, sorry.")
    pass

def uncolor_other_chains():
  mol_id=active_residue()[0]
  ch_id_here=active_residue()[1]  
  blank_colour=0
  blank_res_list=[]
  for ch_id in chain_ids(mol_id):
    sn_max=chain_n_residues(ch_id,mol_id)
    if ch_id!=ch_id_here:
      for sn in range(0,sn_max+1):
        resn=seqnum_from_serial_number(mol_id,ch_id,sn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
        residue_to_color=[([ch_id,resn,ins_id],blank_colour)]
        blank_res_list=blank_res_list+residue_to_color
  try:
    set_user_defined_atom_colour_by_residue_py(mol_id,blank_res_list)
    graphics_to_user_defined_atom_colours_representation(mol_id)
  except NameError:
    info_dialog("You need a newer Coot - custom coloring is only in r6174 and later, sorry.")
    pass

def color_active_chain():
  mol_id=active_residue()[0]
  ch_id_here=active_residue()[1]  
  blank_colour=0
  chain_colour=22 #yellow
  blank_res_list=[]
  chain_res_list=[]
  for ch_id in chain_ids(mol_id):
    sn_max=chain_n_residues(ch_id,mol_id)
    if ch_id!=ch_id_here:
      for sn in range(0,sn_max+1):
        resn=seqnum_from_serial_number(mol_id,ch_id,sn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
        residue_to_color=[([ch_id,resn,ins_id],blank_colour)]
        blank_res_list=blank_res_list+residue_to_color
    if ch_id==ch_id_here:
      for sn in range(0,sn_max+1):
        resn=seqnum_from_serial_number(mol_id,ch_id,sn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
        residue_to_color=[([ch_id,resn,ins_id],chain_colour)]
        chain_res_list=chain_res_list+residue_to_color
  try:
    set_user_defined_atom_colour_by_residue_py(mol_id,blank_res_list)
    set_user_defined_atom_colour_by_residue_py(mol_id,chain_res_list)
    graphics_to_user_defined_atom_colours_representation(mol_id)
  except NameError:
    info_dialog("You need a newer Coot - custom coloring is only in r6174 and later, sorry.")
    pass

def color_active_chain_by_num(chain_colour):
  mol_id=active_residue()[0]
  ch_id_here=active_residue()[1]
  blank_colour=0
  blank_res_list=[]
  chain_res_list=[]
  for ch_id in chain_ids(mol_id):
    sn_max=chain_n_residues(ch_id,mol_id)
    if ch_id!=ch_id_here:
      for sn in range(0,sn_max+1):
        resn=seqnum_from_serial_number(mol_id,ch_id,sn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
        residue_to_color=[([ch_id,resn,ins_id],blank_colour)]
        blank_res_list=blank_res_list+residue_to_color
    if ch_id==ch_id_here:
      for sn in range(0,sn_max+1):
        resn=seqnum_from_serial_number(mol_id,ch_id,sn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
        residue_to_color=[([ch_id,resn,ins_id],chain_colour)]
        chain_res_list=chain_res_list+residue_to_color
  try:
    set_user_defined_atom_colour_by_residue_py(mol_id,blank_res_list)
    set_user_defined_atom_colour_by_residue_py(mol_id,chain_res_list)
    graphics_to_user_defined_atom_colours_representation(mol_id)
  except NameError:
    info_dialog("You need a newer Coot - custom coloring is only in r6174 and later, sorry.")
    pass
    
def color_protein_na(mol_id):
  blank_colour=0
  protein_colour=22 #yellow
  na_colour=31
  protein_res_list=[]
  na_res_list=[]
  blank_res_list=[]
  for ch_id in chain_ids(mol_id):
    sn_max=chain_n_residues(ch_id,mol_id)
    if is_nucleotide_chain_p(mol_id,ch_id):
      for sn in range(0,sn_max+1):
        resn=seqnum_from_serial_number(mol_id,ch_id,sn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
        residue_to_color=[([ch_id,resn,ins_id],na_colour)]
        na_res_list=na_res_list+residue_to_color
    if is_protein_chain_p(mol_id,ch_id):
      for sn in range(0,sn_max+1):
        resn=seqnum_from_serial_number(mol_id,ch_id,sn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
        residue_to_color=[([ch_id,resn,ins_id],protein_colour)]
        protein_res_list=protein_res_list+residue_to_color
    else:
      for sn in range(0,sn_max+1):
        resn=seqnum_from_serial_number(mol_id,ch_id,sn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
        residue_to_color=[([ch_id,resn,ins_id],blank_colour)]
        blank_res_list=blank_res_list+residue_to_color
  try:
    set_user_defined_atom_colour_by_residue_py(mol_id,blank_res_list)
    set_user_defined_atom_colour_by_residue_py(mol_id,protein_res_list)
    set_user_defined_atom_colour_by_residue_py(mol_id,na_res_list)
    graphics_to_user_defined_atom_colours_representation(mol_id)
  except NameError:
    info_dialog("You need a newer Coot - custom coloring is only in r6174 and later, sorry.")
    pass


def color_waters(mol_id):
  water_colour=31
  blank_colour=0
  water_list=[]
  blank_res_list=[]
  for ch_id in chain_ids(mol_id):
    sn_max=chain_n_residues(ch_id,mol_id)
    for sn in range(0,sn_max+1):
      resname_here=resname_from_serial_number(mol_id,ch_id,sn)
      if resname_here=="HOH":
        resn=seqnum_from_serial_number(mol_id,ch_id,sn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
        residue_to_color=[([ch_id,resn,ins_id],water_colour)]
        water_list=water_list+residue_to_color
      else:
        resn=seqnum_from_serial_number(mol_id,ch_id,sn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
        residue_to_color=[([ch_id,resn,ins_id],blank_colour)]
        blank_res_list=blank_res_list+residue_to_color
  try:
    set_user_defined_atom_colour_by_residue_py(mol_id,water_list)
    set_user_defined_atom_colour_by_residue_py(mol_id,blank_res_list)
    graphics_to_user_defined_atom_colours_all_atoms_representation(mol_id)
  except NameError:
    info_dialog("You need a newer Coot - custom coloring is only in r6174 and later, sorry.")
    pass

  
#Search PDB by active chain
#Modify example below. Need to get (clean - waters etc removed) sequence of 
#active chain, search PDB using REST (look up RCSB REST for details), parse output line by line (try using first four characters of each line, or everything between numeric and :, 
#load PDBs (get_ebi_pdb pdb_id) and SSM superpose (superpose imol1 imol2 move_imol2_flag;
#superpose-with-chain-selection imol1 imol2 chain_imol1 chain_imol2 chain_used_flag_imol1 chain_used_flag_imol2 
#move_imol2_copy_flag) each on active chain. Would probably be useful to optionally only load the best matching chain
#of each struc; but equally might be good to keep everything, so that ligands, conserved waters, binding partners etc can be seen.
# should probably restrict to top 10, and exclude those that are 90%+identical... or maybe take min_idpct and max_idpct as parameters?
#
# def query_pdb_by_active_chain():
#   import urllib2
#   url = 'http://www.rcsb.org/pdb/rest/search'
#   mol_id=active_residue()[0]
#   ch_id=active_residue()[1]
#   seq=print_sequence_chain(mol_id,ch_id) # print_seq won't work, as does not return seq but pipes to stdout. need to write func.
#   print seq
#   seq.replace("X","")
#   print seq
#   
#   queryText = """
# <?xml version="1.0" encoding="UTF-8"?>
# <orgPdbQuery>
# <queryType>org.pdb.query.simple.SequenceQuery</queryType>
# <sequence>{sequence}</sequence>
# <eCutOff>0.00001</eCutOff>
# <searchTool>blast</searchTool>
# <sequenceIdentityCutoff>50</sequenceIdentityCutoff>
# </orgPdbQuery>
#   """.format(sequence=seq)
#   
#   print "query:\n", queryText
#   print "querying PDB...\n"
#   req = urllib2.Request(url, data=queryText)
#   f = urllib2.urlopen(req)
#   result = f.read()
#   if result:
#     print "Found number of PDB entries:", result.count('\n')
#     print result
#   else:
#     print "Failed to retrieve results" 
  
#Highlight various items with ball-and-stick and lines
def highlight_chain_breaks():
  mol_id=active_residue()[0]
  clear_ball_and_stick(mol_id)
  turn_off_backup(mol_id)
  obj_number=generic_object_with_name("chain_breaks_{mol_id}".format(mol_id=mol_id))
  generic_object_clear(obj_number)
  protein_resnames=['ALA','UNK','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','ILE','LEU','LYS','MET','MSE','PHE','PRO','SER','THR','TRP','TYR','VAL']
  na_resnames=['A','C','T','G','U']
  missing_segments_list=[]
  for ch_id in chain_ids(mol_id):
    for resn in range(0,chain_n_residues(ch_id,mol_id)):
      resname=resname_from_serial_number(mol_id,ch_id,resn)
      print(resname)
      if (resname in protein_resnames):
        if ((is_term_type_mc_sn(mol_id,ch_id,resn)==1) or (is_term_type_mn_sn(mol_id,ch_id,resn)==1)):
          seqnum=seqnum_from_serial_number(mol_id,ch_id,resn)
#         sel_string="//"+str(ch_id)+"/"+str(seqnum)+"/CA"
#         make_ball_and_stick(mol_id,sel_string,0,0.5,1
          if (is_term_type_mn_sn(mol_id,ch_id,resn)==1):
            x_mn=atom_specs(mol_id,ch_id,seqnum,"","CA","")[-3]
            y_mn=atom_specs(mol_id,ch_id,seqnum,"","CA","")[-2]
            z_mn=atom_specs(mol_id,ch_id,seqnum,"","CA","")[-1]
            x_mid=(x_mn+x_here)/2
            y_mid=(y_mn+y_here)/2
            z_mid=(z_mn+z_here)/2
            res_missing=seqnum-seqnum_from_serial_number(mol_id,ch_id,resn-1)
            #formula: sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)
            distance=((x_here-x_mn)**2+(y_here-y_mn)**2+(z_here-z_mn)**2)**0.5
            print(("distance",distance))
            distance_per_residue=distance/res_missing
            print(("per residue distance",distance_per_residue))
            label_string="{ch_id}  {res_start}...{res_end} ({res_missing} missing residues)".format(ch_id=ch_id,res_start=seqnum_from_serial_number(mol_id,ch_id,resn-1),res_end=seqnum,res_missing=res_missing)
            if distance_per_residue>3.8:
              label_string="MIND THE GAP! {ch_id}  {res_start}...{res_end} ({res_missing} missing residues for {distance:6.1f} A)".format(ch_id=ch_id,res_start=seqnum_from_serial_number(mol_id,ch_id,resn-1),res_end=seqnum,res_missing=res_missing,distance=distance)
            if res_missing >=50:
              label_string="!!!  {ch_id}  {res_start}...{res_end} ({res_missing} missing residues)".format(ch_id=ch_id,res_start=seqnum_from_serial_number(mol_id,ch_id,resn-1),res_end=seqnum,res_missing=res_missing)
            list_entry=[label_string,x_mid,y_mid,z_mid]
            missing_segments_list.append(list_entry)
            if res_missing <=15:
              gap_color="gray"
            elif (res_missing > 15) and (res_missing<50):
              gap_color="orange"
            else:
              gap_color="red"
            if distance_per_residue>3.8:
              gap_color="cyan" 
            line_width=4
            dash_density=3
            to_generic_object_add_dashed_line(obj_number,gap_color,line_width,dash_density,x_here,y_here,z_here,x_mn,y_mn,z_mn)
            #if res_missing>=20:
            #  place_text(str(res_missing),x_mid,y_mid,z_mid,1)
          else: #By definition we always pass through here first (need to hit a term_type_mc befor term_type_mn)
            x_here=atom_specs(mol_id,ch_id,seqnum,"","CA","")[-3]
            y_here=atom_specs(mol_id,ch_id,seqnum,"","CA","")[-2]
            z_here=atom_specs(mol_id,ch_id,seqnum,"","CA","")[-1]
        if (resn==0):
          seqnum=seqnum_from_serial_number(mol_id,ch_id,resn)
          sel_string="//"+str(ch_id)+"/"+str(seqnum)+"/CA"
          make_ball_and_stick(mol_id,sel_string,0,0.5,1)
        if (is_last_polymer_residue_sn(mol_id,ch_id,resn)==1):
          seqnum=seqnum_from_serial_number(mol_id,ch_id,resn)
          sel_string="//"+str(ch_id)+"/"+str(seqnum)+"/CA"
          make_ball_and_stick(mol_id,sel_string,0,0.5,1)
      if (resname in na_resnames):
        if ((is_term_type_mc_sn(mol_id,ch_id,resn)==1) or (is_term_type_mn_sn(mol_id,ch_id,resn)==1)):
          seqnum=seqnum_from_serial_number(mol_id,ch_id,resn)
#         sel_string="//"+str(ch_id)+"/"+str(seqnum)+"/CA"
#         make_ball_and_stick(mol_id,sel_string,0,0.5,1
          if (is_term_type_mn_sn(mol_id,ch_id,resn)==1):
            x_mn=atom_specs(mol_id,ch_id,seqnum,"","P","")[-3]
            y_mn=atom_specs(mol_id,ch_id,seqnum,"","P","")[-2]
            z_mn=atom_specs(mol_id,ch_id,seqnum,"","P","")[-1]
            x_mid=(x_mn+x_here)/2
            y_mid=(y_mn+y_here)/2
            z_mid=(z_mn+z_here)/2
            res_missing=seqnum-seqnum_from_serial_number(mol_id,ch_id,resn-1)
            #formula: sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)
            distance=((x_here-x_mn)**2+(y_here-y_mn)**2+(z_here-z_mn)**2)**0.5
            distance_per_residue=distance/res_missing
            label_string="{ch_id}  {res_start}...{res_end} ({res_missing} missing residues)".format(ch_id=ch_id,res_start=seqnum_from_serial_number(mol_id,ch_id,resn-1),res_end=seqnum,res_missing=res_missing)
            if distance_per_residue>5.9:
              label_string="MIND THE GAP! {ch_id}  {res_start}...{res_end} ({res_missing} missing residues for {distance:6.1f} A)".format(ch_id=ch_id,res_start=seqnum_from_serial_number(mol_id,ch_id,resn-1),res_end=seqnum,res_missing=res_missing,distance=distance)
            if res_missing >=50:
              label_string="!!!  {ch_id}  {res_start}...{res_end} ({res_missing} missing residues)".format(ch_id=ch_id,res_start=seqnum_from_serial_number(mol_id,ch_id,resn-1),res_end=seqnum,res_missing=res_missing)
            list_entry=[label_string,x_mid,y_mid,z_mid]
            missing_segments_list.append(list_entry)
            if res_missing <=15:
              gap_color="gray"
            elif (res_missing > 15) and (res_missing<50):
              gap_color="orange"
            else:
              gap_color="red"
            if distance_per_residue>5.9:
              gap_color="cyan" 
            line_width=4
            dash_density=3
            to_generic_object_add_dashed_line(obj_number,gap_color,line_width,dash_density,x_here,y_here,z_here,x_mn,y_mn,z_mn)
            #if res_missing>=20:
            #  place_text(str(res_missing),x_mid,y_mid,z_mid,1)
          else: #By definition we always pass through here first (need to hit a term_type_mc befor term_type_mn)
            x_here=atom_specs(mol_id,ch_id,seqnum,"","P","")[-3]
            y_here=atom_specs(mol_id,ch_id,seqnum,"","P","")[-2]
            z_here=atom_specs(mol_id,ch_id,seqnum,"","P","")[-1]
        if (resn==0):
          seqnum=seqnum_from_serial_number(mol_id,ch_id,resn)
          sel_string="//"+str(ch_id)+"/"+str(seqnum)+"/P"
          make_ball_and_stick(mol_id,sel_string,0,0.5,1)
        if (is_last_polymer_residue_sn(mol_id,ch_id,resn)==1):
          seqnum=seqnum_from_serial_number(mol_id,ch_id,resn)
          sel_string="//"+str(ch_id)+"/"+str(seqnum)+"/P"
          make_ball_and_stick(mol_id,sel_string,0,0.5,1)
  set_display_generic_object(obj_number,1)
  try:
    attach_generic_object_to_molecule(obj_number, mol_id)
  except NameError:
    info_dialog("attach_generic_object_to_molecule is only present in Coot r6057 and later, sorry.")
    pass
  turn_on_backup(mol_id)
  interesting_things_gui("Missing segments",missing_segments_list)


def highlight_all_chain_breaks():
  protein_resnames=['ALA','UNK','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','ILE','LEU','LYS','MET','MSE','PHE','PRO','SER','THR','TRP','TYR','VAL']
  na_resnames=['A','C','T','G','U']
  for mol_id in model_molecule_list():
    clear_ball_and_stick(mol_id)
    turn_off_backup(mol_id)
    obj_number=generic_object_with_name("chain_breaks_{mol_id}".format(mol_id=mol_id))
    generic_object_clear(obj_number)
    for ch_id in chain_ids(mol_id):
      for resn in range(0,chain_n_residues(ch_id,mol_id)):
        resname=resname_from_serial_number(mol_id,ch_id,resn)
        print(resname)
        if (resname in protein_resnames):
          if ((is_term_type_mc_sn(mol_id,ch_id,resn)==1) or (is_term_type_mn_sn(mol_id,ch_id,resn)==1)):
            seqnum=seqnum_from_serial_number(mol_id,ch_id,resn)
  #         sel_string="//"+str(ch_id)+"/"+str(seqnum)+"/CA"
  #         make_ball_and_stick(mol_id,sel_string,0,0.5,1)
            if (is_term_type_mn_sn(mol_id,ch_id,resn)==1):
              x_mn=atom_specs(mol_id,ch_id,seqnum,"","CA","")[-3]
              y_mn=atom_specs(mol_id,ch_id,seqnum,"","CA","")[-2]
              z_mn=atom_specs(mol_id,ch_id,seqnum,"","CA","")[-1]
              x_mid=(x_mn+x_here)/2
              y_mid=(y_mn+y_here)/2
              z_mid=(z_mn+z_here)/2
              res_missing=seqnum-seqnum_from_serial_number(mol_id,ch_id,resn-1)
              #formula: sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)
              distance=((x_here-x_mn)**2+(y_here-y_mn)**2+(z_here-z_mn)**2)**0.5
              distance_per_residue=distance/res_missing
              if res_missing <=15:
                gap_color="gray"
              elif (res_missing > 15) and (res_missing<50):
                gap_color="orange"
              else:
                gap_color="red"
              if distance_per_residue>3.8:
                gap_color="cyan" 
              line_width=4
              dash_density=3
              to_generic_object_add_dashed_line(obj_number,gap_color,line_width,dash_density,x_here,y_here,z_here,x_mn,y_mn,z_mn)
              #res_missing=seqnum-seqnum_from_serial_number(mol_id,ch_id,resn-1)
              #if res_missing>=20:
              #  place_text(str(res_missing),x_mid,y_mid,z_mid,1)
            else: #By definition we always pass through here first (need to hit a term_type_mc befor term_type_mn)
              x_here=atom_specs(mol_id,ch_id,seqnum,"","CA","")[-3]
              y_here=atom_specs(mol_id,ch_id,seqnum,"","CA","")[-2]
              z_here=atom_specs(mol_id,ch_id,seqnum,"","CA","")[-1]
          if (resn==0):
            seqnum=seqnum_from_serial_number(mol_id,ch_id,resn)
            sel_string="//"+str(ch_id)+"/"+str(seqnum)+"/CA"
            make_ball_and_stick(mol_id,sel_string,0,0.5,1)
          if (is_last_polymer_residue_sn(mol_id,ch_id,resn)==1):
            seqnum=seqnum_from_serial_number(mol_id,ch_id,resn)
            sel_string="//"+str(ch_id)+"/"+str(seqnum)+"/CA"
            make_ball_and_stick(mol_id,sel_string,0,0.5,1)
        if (resname in na_resnames):
          if ((is_term_type_mc_sn(mol_id,ch_id,resn)==1) or (is_term_type_mn_sn(mol_id,ch_id,resn)==1)):
            seqnum=seqnum_from_serial_number(mol_id,ch_id,resn)
            if (is_term_type_mn_sn(mol_id,ch_id,resn)==1):
              x_mn=atom_specs(mol_id,ch_id,seqnum,"","P","")[-3]
              y_mn=atom_specs(mol_id,ch_id,seqnum,"","P","")[-2]
              z_mn=atom_specs(mol_id,ch_id,seqnum,"","P","")[-1]
              x_mid=(x_mn+x_here)/2
              y_mid=(y_mn+y_here)/2
              z_mid=(z_mn+z_here)/2
              res_missing=seqnum-seqnum_from_serial_number(mol_id,ch_id,resn-1)
              #formula: sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)
              distance=((x_here-x_mn)**2+(y_here-y_mn)**2+(z_here-z_mn)**2)**0.5
              distance_per_residue=distance/res_missing
              if res_missing <=15:
                gap_color="gray"
              elif (res_missing > 15) and (res_missing<50):
                gap_color="orange"
              else:
                gap_color="red"
              if distance_per_residue>5.9:
                gap_color="cyan" 
              line_width=4
              dash_density=3
              to_generic_object_add_dashed_line(obj_number,gap_color,line_width,dash_density,x_here,y_here,z_here,x_mn,y_mn,z_mn)
              #if res_missing>=20:
              #  place_text(str(res_missing),x_mid,y_mid,z_mid,1)
            else: #By definition we always pass through here first (need to hit a term_type_mc befor term_type_mn)
              x_here=atom_specs(mol_id,ch_id,seqnum,"","P","")[-3]
              y_here=atom_specs(mol_id,ch_id,seqnum,"","P","")[-2]
              z_here=atom_specs(mol_id,ch_id,seqnum,"","P","")[-1]
          if (resn==0):
            seqnum=seqnum_from_serial_number(mol_id,ch_id,resn)
            sel_string="//"+str(ch_id)+"/"+str(seqnum)+"/P"
            make_ball_and_stick(mol_id,sel_string,0,0.5,1)
          if (is_last_polymer_residue_sn(mol_id,ch_id,resn)==1):
            seqnum=seqnum_from_serial_number(mol_id,ch_id,resn)
            sel_string="//"+str(ch_id)+"/"+str(seqnum)+"/P"
            make_ball_and_stick(mol_id,sel_string,0,0.5,1)
    set_display_generic_object(obj_number,1) 
    try:
      attach_generic_object_to_molecule(obj_number, mol_id)
    except NameError:
      info_dialog("attach_generic_object_to_molecule is only present in Coot r6057 and later, sorry.")
      pass
    turn_on_backup(mol_id)

def local_hbonds_and_baddies():
  mol_id=active_residue()[0]
  probe_local_sphere(mol_id,5)
  undisplay_list=["wide contact","close contact"]
  display_list=["bad overlap","H-bonds","small overlap"]
  for obj_name in undisplay_list:
    obj_number=generic_object_with_name(obj_name)
    generic_object_clear(obj_number)
  for obj_name in display_list:
    obj_number=generic_object_with_name(obj_name)
    set_display_generic_object(obj_number,1)
    try:
      attach_generic_object_to_molecule(obj_number, mol_id)
    except NameError:
      info_dialog("attach_generic_object_to_molecule is only present in Coot r6057 and later, sorry.")
      pass
      
def global_hbonds_and_baddies():
  mol_id=active_residue()[0]
  probe(mol_id)
  close_molecule(model_molecule_list()[-1])
  undisplay_list=["wide contact","close contact","H-bonds","small overlap"]
  display_list=["bad overlap"]
  for obj_name in undisplay_list:
    obj_number=generic_object_with_name(obj_name)
    generic_object_clear(obj_number)
  for obj_name in display_list:
    obj_number=generic_object_with_name(obj_name)
    set_display_generic_object(obj_number,1)
    try:
      attach_generic_object_to_molecule(obj_number, mol_id)
    except NameError:
      info_dialog("attach_generic_object_to_molecule is only present in Coot r6057 and later, sorry.")
      pass
      
def clear_dots():
  undisplay_list=["wide contact","close contact", "small overlap","bad overlap","H-bonds","H-bond","big-overlap","close-contact","small-overlap","wide-contact","clashes"]
  for obj_name in undisplay_list:
    obj_number=generic_object_with_name(obj_name)
    generic_object_clear(obj_number)

#Place helix and add helix restraints
def place_helix_with_restraints():
  place_helix_here()

#Merge two chains (throw error msg if they overlap in numbering)
#This needs fixing - says selections overlap in cases where sel1 has one segment on either side from sel2
#Fixed now, I think.
def merge_chains():
  def merge_chains_by_click(res1,res2):
    mol_id_1=res1[1]
    mol_id_2=res2[1]
    ch_id_1=res1[2]
    ch_id_2=res2[2]
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
  else:
    mol_id=active_residue()[0]
    segments=segment_list(mol_id)
    for seg in segments:
      res_start=seg[2]
      res_end=seg[3]
      ch_id=seg[1]
      turn_off_backup(mol_id)
      set_refinement_immediate_replacement(1)
      rigid_body_refine_zone(res_start,res_end,ch_id,mol_id)
      accept_regularizement()
      set_refinement_immediate_replacement(0)
      turn_on_backup(mol_id)
      
#Fit current segment
def fit_this_segment():
  if (imol_refinement_map()==-1):
    info_dialog("You must set a refinement map!")
  else: 
    mol_id=active_residue()[0]
    segments=segment_list(mol_id)
    res_here=active_residue()[2]
    ch_id=active_residue()[1]
    for seg in segments:
      if (res_here>=seg[2]) and (res_here<=seg[3]) and (ch_id==seg[1]):
        res_start=seg[2]
        res_end=seg[3]
        ch_id=seg[1]
        turn_off_backup(mol_id)
        set_refinement_immediate_replacement(1)
        rigid_body_refine_zone(res_start,res_end,ch_id,mol_id)
        accept_regularizement()
        set_refinement_immediate_replacement(0)
        turn_on_backup(mol_id)

#Set default b-fac for new atoms to mean B for active mol
def set_new_atom_b_fac_to_mean():
  mol_id=active_residue()[0]
  mean_b=average_temperature_factor(mol_id)
  set_default_temperature_factor_for_new_atoms(mean_b)
  
  
#Shortcut for adding terminal residue (to bind to key)
def add_term_shortcut():
  if imol_refinement_map()!=-1:
    mol_id=active_residue()[0]
    ch_id=active_residue()[1]
    resn=active_residue()[2]
    atom_name=active_residue()[4]
    first_in_seg=first_residue_in_seg(mol_id,ch_id,resn)
    last_in_seg=last_residue_in_seg(mol_id,ch_id,resn)
    delta_first=abs(first_in_seg-resn)
    delta_last=abs(last_in_seg-resn)
    set_new_atom_b_fac_to_mean()
    if delta_first<=delta_last:
      set_go_to_atom_molecule(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      add_terminal_residue(mol_id,ch_id,first_in_seg,"auto",1)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg-1,atom_name)
    else:
      set_go_to_atom_molecule(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      add_terminal_residue(mol_id,ch_id,last_in_seg,"auto",1)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg+1,atom_name)
  else:
    info_dialog("You must set a refinement map!")
    
def add_term_shortcut_force():
  mol_id=active_residue()[0]
  ch_id=active_residue()[1]
  resn=active_residue()[2]
  atom_name=active_residue()[4]
  set_go_to_atom_molecule(mol_id)
  first_in_seg=first_residue_in_seg(mol_id,ch_id,resn)
  last_in_seg=last_residue_in_seg(mol_id,ch_id,resn)
  delta_first=abs(first_in_seg-resn)
  delta_last=abs(last_in_seg-resn)
  set_new_atom_b_fac_to_mean()
  if delta_first<=delta_last:
    set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
    force_add_terminal_residue_noclick(mol_id,ch_id,first_in_seg)
    sort_residues(mol_id)
    set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg-1,atom_name)
  else:
    set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
    force_add_terminal_residue_noclick(mol_id,ch_id,last_in_seg)
    sort_residues(mol_id)
    set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg+1,atom_name)
    
def add_term_shortcut_force_strand():
  mol_id=active_residue()[0]
  ch_id=active_residue()[1]
  resn=active_residue()[2]
  atom_name=active_residue()[4]
  set_go_to_atom_molecule(mol_id)
  first_in_seg=first_residue_in_seg(mol_id,ch_id,resn)
  last_in_seg=last_residue_in_seg(mol_id,ch_id,resn)
  delta_first=abs(first_in_seg-resn) 
  delta_last=abs(last_in_seg-resn)
  set_new_atom_b_fac_to_mean()
  if delta_first<=delta_last:
    set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
    force_add_terminal_residue_noclick_strand(mol_id,ch_id,first_in_seg)
    sort_residues(mol_id)
    set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg-1,atom_name)
  else:
    set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
    force_add_terminal_residue_noclick_strand(mol_id,ch_id,last_in_seg)
    sort_residues(mol_id)
    set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg+1,atom_name)

#Add h-bond restraints to active mol with Prosmart
def run_prosmart_self():
  """
  target is my molecule, ref is the homologous (high-res) model

  extra arg: include_side_chains=False
  """
  imol_target=active_residue()[0]
  dir_stub = "coot-ccp4"
  prosmart_dir="ProSMART_Output"
  make_directory_maybe(dir_stub)
  make_directory_maybe(prosmart_dir)
  target_pdb_file_name = os.path.join(dir_stub,
                                      molecule_name_stub(imol_target, 0).replace(" ", "_") + \
                                      "-prosmart.pdb")
  prosmart_out = os.path.join("ProSMART_Output",
                              molecule_name_stub(imol_target, 0).replace(" ", "_") + \
                              "-prosmart.txt")

  write_pdb_file(imol_target, target_pdb_file_name)
  prosmart_exe = find_exe("prosmart")
  if prosmart_exe:
    l = ["-p1", target_pdb_file_name,
         "-h","-self_restrain","-bond_max","4.0","-bond_override", "2","-o",prosmart_dir]
    popen_command(prosmart_exe,
                  l,
                  [],
                  os.path.join(dir_stub, "prosmart.log"),
                  False)
    if (not os.path.isfile(prosmart_out)):
      print("file not found", prosmart_out)
    else:
      print("Reading ProSMART restraints from", prosmart_out)
      add_refmac_extra_restraints(imol_target, prosmart_out)
      delete_extra_restraints_worse_than(imol_target,4)
      
#Fix all atoms in active residue
def fix_active_residue():
  ar=active_residue() #ar[0] is mol_id, 1 is ch_id, 2 is resnum, 3 is ins code, 4 is atom name (CA) and 5 is alt conf 
  residue_info_list=residue_info(ar[0],ar[1],ar[2],ar[3])
  atom_list=[]
  for item in residue_info_list: # residue item is list of list of lists. first item in first list of each list in the list is the atom name :-)
    atom_name=item[0][0]
    atom_info=[ar[1],ar[2],ar[3],atom_name,ar[5]] #Make a list of chain id, res number, ins code, atom name and alt conf. Not sure if alt conf and ins code are right way around here...
    atom_list.append(atom_info)
  mark_multiple_atoms_as_fixed(ar[0],atom_list,1)
  
#Save and overwrite active model:
def quicksave_active():
  import shutil
  mol_id=active_residue()[0]
  filename=molecule_name(mol_id)
  filename_bak=filename+".bak"
  shutil.copy(filename,filename_bak)
  save_coordinates(mol_id,filename)
  info_dialog("Saved {filename} (molecule {mol_id}). \n If this was an accident, you can find a backup of the original file here: \n {filename_bak}".format(filename=filename,mol_id=mol_id,filename_bak=filename_bak))

#Refresh model from file (if changed externally for example)
def reload_model():
  mol_id=active_residue()[0]
  filename=molecule_name(mol_id)
  clear_and_update_model_molecule_from_file(mol_id,filename)

#Copy changes from active chain to NCS equivalents.
def copy_ncs_chain_from_active():
  mol_id=active_residue()[0]
  ch_id=active_residue()[1]
  turn_off_backup(mol_id)
  ncs_control_change_ncs_master_to_chain_id(mol_id,ch_id)
  copy_from_ncs_master_to_others(mol_id,ch_id)
  turn_on_backup(mol_id)

def return_seq_as_string(mol_id,ch_id):
  seq=''
  for sn in range(0,chain_n_residues(ch_id,mol_id)):
    resno=seqnum_from_serial_number(mol_id,ch_id,sn)
    ins_code=insertion_code_from_serial_number(mol_id,ch_id,sn)
    residue_name_3_letter=residue_name(mol_id,ch_id,resno,ins_code)
    aa_code=three_letter_code2single_letter(residue_name_3_letter)
    if (len(residue_name(mol_id,ch_id,resno,ins_code))==1):
      aa_code=residue_name_3_letter
    if (residue_name(mol_id,ch_id,resno,ins_code)!="ALA") and (aa_code=="A") and (len(residue_name(mol_id,ch_id,resno,ins_code))!=1):
      aa_code="X"
    seq=seq+aa_code
  return seq

stereo_gif_counter=1
def make_rot_gif():
  if find_exe("convert"):
    global stereo_gif_counter
    import subprocess
    pwd=os.getcwd()
    set_rotation_center_size(0)
    set_draw_axes(0)
    rotate_y_scene(1,-2)
    screendump_image("{pwd}/pair1_tmp.ppm".format(pwd=pwd))
    rotate_y_scene(1,1)
    screendump_image("{pwd}/pair2_tmp.ppm".format(pwd=pwd))
    rotate_y_scene(1,1)
    screendump_image("{pwd}/pair3_tmp.ppm".format(pwd=pwd))
    rotate_y_scene(1,1)
    screendump_image("{pwd}/pair4_tmp.ppm".format(pwd=pwd))
    rotate_y_scene(1,1)
    screendump_image("{pwd}/pair5_tmp.ppm".format(pwd=pwd))
    rotate_y_scene(1,-2)
    p=subprocess.Popen("convert -scale 500 -normalize -fuzz 5% -delay 8 -loop 0 -layers Optimize {pwd}/pair1_tmp.ppm {pwd}/pair2_tmp.ppm {pwd}/pair3_tmp.ppm {pwd}/pair4_tmp.ppm {pwd}/pair5_tmp.ppm {pwd}/pair4_tmp.ppm {pwd}/pair3_tmp.ppm {pwd}/pair2_tmp.ppm {pwd}/stereo_density_{stereo_gif_counter}.gif".format(pwd=pwd,stereo_gif_counter=stereo_gif_counter),shell=True) 
    p.communicate()
    stereo_gif_counter=stereo_gif_counter+1
    os.remove("{pwd}/pair1_tmp.ppm".format(pwd=pwd))
    os.remove("{pwd}/pair2_tmp.ppm".format(pwd=pwd))
    os.remove("{pwd}/pair3_tmp.ppm".format(pwd=pwd))
    os.remove("{pwd}/pair4_tmp.ppm".format(pwd=pwd))
    os.remove("{pwd}/pair5_tmp.ppm".format(pwd=pwd))
    set_rotation_center_size(0.1)
    set_draw_axes(1)
  else:
    info_dialog("You need Imagemagick to use this!")
    
    
def find_sequence_in_current_chain(subseq):
  subseq=subseq.upper()
  mol_id=active_residue()[0]
  ch_id=active_residue()[1]
  seq=return_seq_as_string(mol_id,ch_id)
  index=0
  sn_list=[]
  interesting_list=[]
  while index < len(seq):
    try:
      index=seq.index(subseq,index)
    except ValueError:
      break
    sn_list.append(index)
    if index==-1:
      break
    index=index+len(subseq)
  if len(sn_list)==1:
    sn_start=sn_list[0]
    resno=seqnum_from_serial_number(mol_id,ch_id,sn_start)
    ins_code=insertion_code_from_serial_number(mol_id,ch_id,sn_start)
    alt_conf=residue_alt_confs(mol_id,ch_id,resno,ins_code)[0]
    print(("sn_start",sn_start))
    try:
      x=atom_specs(mol_id,ch_id,resno,ins_code,"CA",alt_conf)[-3]
      y=atom_specs(mol_id,ch_id,resno,ins_code,"CA",alt_conf)[-2]
      z=atom_specs(mol_id,ch_id,resno,ins_code,"CA",alt_conf)[-1]
    except TypeError:
      x=atom_specs(mol_id,ch_id,resno,ins_code,"P",alt_conf)[-3]
      y=atom_specs(mol_id,ch_id,resno,ins_code,"P",alt_conf)[-2]
      z=atom_specs(mol_id,ch_id,resno,ins_code,"P",alt_conf)[-1]
    set_rotation_centre(x,y,z)
  elif len(sn_list)==0:
    info_dialog("Sequence not found!")
  elif len(sn_list)>1:
    count=0
    for sn in sn_list:
      count=count+1
      resno=seqnum_from_serial_number(mol_id,ch_id,sn)
      ins_code=insertion_code_from_serial_number(mol_id,ch_id,sn)
      alt_conf=residue_alt_confs(mol_id,ch_id,resno,ins_code)[0]
      print(("sn_start",sn))
      try:
        x=atom_specs(mol_id,ch_id,resno,ins_code,"CA",alt_conf)[-3]
        y=atom_specs(mol_id,ch_id,resno,ins_code,"CA",alt_conf)[-2]
        z=atom_specs(mol_id,ch_id,resno,ins_code,"CA",alt_conf)[-1]
      except TypeError:
        x=atom_specs(mol_id,ch_id,resno,ins_code,"P",alt_conf)[-3]
        y=atom_specs(mol_id,ch_id,resno,ins_code,"P",alt_conf)[-2]
        z=atom_specs(mol_id,ch_id,resno,ins_code,"P",alt_conf)[-1]
      list_entry=[str(resno),x,y,z]
      interesting_list.append(list_entry)
      print(("interesting list",interesting_list))
    print(("interesting list",interesting_list))
    interesting_things_gui("Matches to entered sequence",interesting_list)

#Test post manipulation background func
# import time, threading
# def background_func():
#   threading.Timer(0.01, post_manip_background).start()
#   if condition:
#     do something
# post_manip_background()

# def post_manipulation_script(imol, mode):
#   print "BL DEBUG:: imol and mode", imol, mode
#   if (mode == DELETED):
#     print "BL DEBUG:: deleted something in mol ", imol
#     return 1
# 
# def post_manipulation_script2(imol, mode):
#   print "BL DEBUG:: imol and mode", imol, mode
#   if (mode == MUTATED):
#     print "BL DEBUG:: moved something in mol ", imol
#     return 1

#****Disabled Gtk-dependent UI archive****
# Re-enable these once embedded Python can import Gtk/coot_gui cleanly again.
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
# Disabled Gtk/dialog helper functions are preserved verbatim here so they can
# be restored cleanly once embedded Python can import Gtk/coot_gui again.
GTK_DIALOG_FUNCTION_ARCHIVE = r'''
def set_map_level_quickly():
  if scroll_wheel_map()!=-1 and map_is_displayed(scroll_wheel_map())!=0:
    current_map_level=get_contour_level_in_sigma(scroll_wheel_map())
    current_map_level="{0:.2f}".format(current_map_level)
    def set_map_level_quickly(X):
      try:
        map_level=float(X)
        map_id=scroll_wheel_map()
        set_contour_level_in_sigma(map_id, map_level)
      except ValueError:
        info_dialog("Has to be a number!") 
    generic_single_entry("New map level in sigma/RMS?",current_map_level,"Set map level",set_map_level_quickly)
  else:
    info_dialog("You need a (scrollable, displayed) map!")

#Grow helix from selected terminus
def grow_helix():
  def grow_helix_post_click(res1):
    def grow_helix_enter_resn(n):
      mol_id=res1[1]
      ch_id=res1[2]
      res_no=res1[3]
      res_no_0=res_no
      for i in range(1,(int(n)+1)):
        res_type="auto"
        add_terminal_residue_using_phi_psi(mol_id,ch_id,res_no,
        res_type,-57.82,-47)
        sort_residues(mol_id)
        if (res_no==(first_residue_in_seg(mol_id,ch_id,res_no)+1)):
          res_no=res_no-1
        elif (res_no==(last_residue_in_seg(mol_id,ch_id,res_no)-1)):
          res_no=res_no+1
      set_b_factor_residue_range(mol_id,ch_id,res_no_0,res_no,default_new_atoms_b_factor())
    generic_single_entry("How many residues for helix?",
    "10","Grow helix",grow_helix_enter_resn)
  user_defined_click(1,grow_helix_post_click)
  
#Grow strand from selected terminus
def grow_strand():
  def grow_strand_post_click(res1):
    def grow_strand_enter_resn(n):
      mol_id=res1[1]
      ch_id=res1[2]
      res_no=res1[3]
      res_no_0=res_no
      for i in range(1,(int(n)+1)):
        res_type="auto"
        add_terminal_residue_using_phi_psi(mol_id,ch_id,res_no,
        res_type,-139,135)
        sort_residues(mol_id)
        if (res_no==(first_residue_in_seg(mol_id,ch_id,res_no)+1)):
          res_no=res_no-1
        elif (res_no==(last_residue_in_seg(mol_id,ch_id,res_no)-1)):
          res_no=res_no+1
      set_b_factor_residue_range(mol_id,ch_id,res_no_0,res_no,default_new_atoms_b_factor())
    generic_single_entry("How many residues for strand?",
    "10","Grow strand",grow_strand_enter_resn)
  user_defined_click(1,grow_strand_post_click)
  
#Grow para strand from selected terminus
def grow_parallel_strand():
  def grow_parallel_strand_post_click(res1):
    def grow_parallel_strand_enter_resn(n):
      mol_id=res1[1]
      ch_id=res1[2]
      res_no=res1[3]
      for i in range(1,(int(n)+1)):
        res_type="auto"
        add_terminal_residue_using_phi_psi(mol_id,ch_id,res_no,
        res_type,-119,113)
        sort_residues(mol_id)
        if (res_no==(first_residue_in_seg(mol_id,ch_id,res_no)+1)):
          res_no=res_no-1
        elif (res_no==(last_residue_in_seg(mol_id,ch_id,res_no)-1)):
          res_no=res_no+1
    generic_single_entry("How many residues for parallel strand?",
    "10","Grow parallel strand",grow_parallel_strand_enter_resn)
  user_defined_click(1,grow_parallel_strand_post_click)

#Grow 3-10 helix from selected terminus
def grow_helix_3_10():
  def grow_helix_post_click(res1):
    def grow_helix_enter_resn(n):
      mol_id=res1[1]
      ch_id=res1[2]
      res_no=res1[3]
      for i in range(1,(int(n)+1)):
        res_type="auto"
        add_terminal_residue_using_phi_psi(mol_id,ch_id,res_no,
        res_type,-49,-26)
        sort_residues(mol_id)
        if (res_no==(first_residue_in_seg(mol_id,ch_id,res_no)+1)):
          res_no=res_no-1
        elif (res_no==(last_residue_in_seg(mol_id,ch_id,res_no)-1)):
          res_no=res_no+1
    generic_single_entry("How many residues for helix?",
    "10","Grow 3-10 helix",grow_helix_enter_resn)
  user_defined_click(1,grow_helix_post_click)

#Renumbers the active chain (from active_residue()). User enters new starting residue number.
def renumber_by_first_res():
  def renum_chain(new_num):
    new_num=int(new_num)
    mol_id=active_residue()[0]
    ch_id=active_residue()[1]
    first_res=first_residue(mol_id,ch_id)
    last_res=last_residue(mol_id,ch_id)
    offset=new_num-first_res
    renumber_residue_range(mol_id,ch_id,first_res,last_res,int(offset))
    delete_all_extra_restraints(mol_id)
    set_show_extra_restraints(mol_id,0)
    set_show_extra_restraints(mol_id,1)
  generic_single_entry("New number for first residue in this chain?",
  str(seqnum_from_serial_number(active_residue()[0],
  "%s"%(active_residue()[1]),0)),"Renumber",renum_chain)
  
#Renumbers the active chain (from active_residue()). User enters new last residue number.
def renumber_by_last_res():
  def renum_chain(new_num):
    new_num=int(new_num)
    mol_id=active_residue()[0]
    ch_id=active_residue()[1]
    first_res=first_residue(mol_id,ch_id)
    last_res=last_residue(mol_id,ch_id)
    offset=new_num-last_res
    renumber_residue_range(mol_id,ch_id,first_res,last_res,int(offset))
    delete_all_extra_restraints(mol_id)
    set_show_extra_restraints(mol_id,0)
    set_show_extra_restraints(mol_id,1)
  generic_single_entry("New number for last residue in this chain?",
  str(seqnum_from_serial_number(active_residue()[0],
  "%s"%(active_residue()[1]),
  (chain_n_residues(active_residue()[1],active_residue()[0])-1))),"Renumber",renum_chain)
  
#Renumbers the active chain (from active_residue()). User enters new active residue number.
def renumber_by_active_res():
  def renum_chain(new_num):
    new_num=int(new_num)
    mol_id=active_residue()[0]
    ch_id=active_residue()[1]
    current_num=active_residue()[2]
    first_res=first_residue(mol_id,ch_id)
    last_res=last_residue(mol_id,ch_id)
    offset=new_num-current_num
    renumber_residue_range(mol_id,ch_id,first_res,last_res,int(offset))
    delete_all_extra_restraints(mol_id)
    set_show_extra_restraints(mol_id,0)
    set_show_extra_restraints(mol_id,1)
  generic_single_entry("New number for this residue?",
  str(active_residue()[2]),"Renumber",renum_chain)

#Renumber from N-term to current residue
def renumber_n_term_segment():
  def renumber_n_term_segment_entry(new_resn):
    mol_id=active_residue()[0]
    ch_id=active_residue()[1]
    resn=active_residue()[2]
    first_res=seqnum_from_serial_number(mol_id,"%s"%(ch_id),0)
    offset=int(new_resn)-resn
    renumber_residue_range(mol_id,ch_id,first_res,resn,offset)
    delete_all_extra_restraints(mol_id)
    set_show_extra_restraints(mol_id,0)
    set_show_extra_restraints(mol_id,1)
  generic_single_entry("New residue number?",
  str(active_residue()[2]),"Renumber",renumber_n_term_segment_entry)
  
#Renumber from current residue to C-term
def renumber_c_term_segment():
  def renumber_c_term_segment_entry(new_resn):
    mol_id=active_residue()[0]
    ch_id=active_residue()[1]
    resn=active_residue()[2]
    def last_residue(mol_id,ch_id):
      n=chain_n_residues(ch_id,mol_id)-1
      result=seqnum_from_serial_number(mol_id,"%s"%(ch_id),n)
      return result
    last_res=last_residue(mol_id,ch_id)
    offset=int(new_resn)-resn
    renumber_residue_range(mol_id,ch_id,resn,last_res,offset)
    delete_all_extra_restraints(mol_id)
    set_show_extra_restraints(mol_id,0)
    set_show_extra_restraints(mol_id,1)
  generic_single_entry("New residue number?",
  str(active_residue()[2]),"Renumber",renumber_c_term_segment_entry)
  
#Sharpen or blur map. Slider jerky on large maps.
def sharpen_by_entered_factor():
  if (scroll_wheel_map()==-1):
    info_dialog("You need a map!")
  else:
    def sharpen_by_entry(B):
      B=float(B)
      sharpen(scroll_wheel_map(),B)
    generic_single_entry("New B-factor for map?",
    "-100","Sharpen/Blur",sharpen_by_entry)
    
#Reload map with different hi-res limit
def change_hires_limit():
  if (scroll_wheel_map()==-1):
    info_dialog("You need a map!")
  else:
    def change_hires_by_entry(new_res):
      mol=scroll_wheel_map()
      new_res=float(new_res)
      mtz_file=map_parameters(mol)[0]
      F_col=map_parameters(mol)[1]
      PHI_col=map_parameters(mol)[2]
      make_and_draw_map_with_reso_with_refmac_params(mtz_file,F_col,PHI_col,"",0,0,0,"Fobs:None-specified",
      "SigF:None-specified","RFree:None-specified",0,0,1,1000.0,new_res)
      close_molecule(mol)
    generic_single_entry("New high-res limit for map?",
    "5.0","Change high resolution limit for active map",change_hires_by_entry)
    
def change_hires_limit_copy():
  if (scroll_wheel_map()==-1):
    info_dialog("You need a map!")
  else:
    def change_hires_by_entry(new_res):
      mol=scroll_wheel_map()
      new_res=float(new_res)
      mtz_file=map_parameters(mol)[0]
      F_col=map_parameters(mol)[1]
      PHI_col=map_parameters(mol)[2]
      make_and_draw_map_with_reso_with_refmac_params(mtz_file,F_col,PHI_col,"",0,0,0,"Fobs:None-specified",
      "SigF:None-specified","RFree:None-specified",0,0,1,1000.0,new_res)
    generic_single_entry("New high-res limit for map?",
    "5.0","Make low pass filtered copy of active map",change_hires_by_entry)

#Make a button list for inserting common monomers
def pick_common_monomers():
  get_ddm=["DDM", lambda func: get_monomer_no_H("LMT")]
  get_dm=["DM", lambda func: get_monomer_no_H("DMU")]
  get_bog=["OG", lambda func: get_monomer_no_H("HSH")]
  get_ldao=["LDAO", lambda func: get_monomer_no_H("LDA")]
  get_mpg=["Monoolein", lambda func: get_monomer_no_H("MPG")]
  get_glycerol=["Glycerol", lambda func: get_monomer_no_H("GOL")]
  get_eg=["Ethylene glycol", lambda func: get_monomer_no_H("EDO")]
  get_acetate=["Acetate", lambda func: get_monomer_no_H("ACT")]
  get_dmso=["DMSO", lambda func: get_monomer_no_H("DMS")]
  get_tris=["Tris", lambda func: get_monomer_no_H("TAM")]
  get_hepes=["HEPES", lambda func: get_monomer_no_H("EPE")]
  get_mes=["MES", lambda func: get_monomer_no_H("MES")]
  get_cac=["Cacodylate", lambda func: get_monomer_no_H("CAD")]
  get_peg=["PEG", lambda func: get_monomer_no_H("1PE")]
  get_popg=["POPG (lipid)", lambda func: get_monomer_no_H("LHG")]
  get_pe=["PE (lipid)", lambda func: get_monomer_no_H("PEF")]
  get_pc=["PC (lipid)", lambda func: get_monomer_no_H("PLC")]
  get_ps=["PS (lipid)", lambda func: get_monomer_no_H("PSF")]
  get_ps=["PI(3,4)P2 (lipid)", lambda func: get_monomer_no_H("52N")]
  get_chs=["Cholesterol hemisuccinate", lambda func: get_monomer_no_H("Y01")]
  get_chl=["Cholesterol", lambda func: get_monomer_no_H("CLR")]
  button_list=[get_acetate,get_eg,get_glycerol,get_dmso,get_ddm,get_dm,get_bog,get_ldao,get_mpg,get_tris,get_hepes,get_mes,get_cac,get_peg,get_popg,get_pe,get_pc,get_ps,get_chs,get_chl]
  generic_button_dialog("Common small molecules",button_list)

#Colors subset of protein residues red, provided by user as string of single-letter ids.
def color_protein_residue_subset():
  def color_from_string(X):
    entry=str(X).upper() #capitalize
    mol_id=active_residue()[0]
    aa_dic={'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS','E':'GLU','Q':'GLN','G':'GLY','H':'HIS','I':'ILE','L':'LEU','K':'LYS','M':'MET','F':'PHE','P':'PRO','S':'SER','T':'THR','W':'TRP','Y':'TYR','V':'VAL','X':'UNK'}
    highlight_colour=34 #red
    blank_colour=0 #light grey
    residue_list=[]
    resname_list=[]
    for char in entry:
      if (aa_dic.get(char,0)!=0): #make a list of 3-letter resnames from user supplied string by checking aa_dic
        resname=[(aa_dic.get(char,0))]
        resname_list=resname_list+resname
    for ch_id in chain_ids(mol_id):
      sn_max=chain_n_residues(ch_id,mol_id)
      for sn in  range(0,sn_max+1):
        resname_here=resname_from_serial_number(mol_id,ch_id,sn)
        for resname in resname_list:
          if resname==resname_here:
            resn=seqnum_from_serial_number(mol_id,ch_id,sn)
            ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
            residue_to_color=[([ch_id,resn,ins_id],highlight_colour)]
            residue_list=residue_list+residue_to_color
        if resname_here not in resname_list:
          resn=seqnum_from_serial_number(mol_id,ch_id,sn)
          ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
          residue_to_color=[([ch_id,resn,ins_id],blank_colour)]
          residue_list=residue_list+residue_to_color
    try:
      set_user_defined_atom_colour_by_residue_py(mol_id,residue_list)
      graphics_to_user_defined_atom_colours_representation(mol_id)
    except NameError:
      info_dialog("You need a newer Coot - custom coloring is only in r6174 and later, sorry.")
      pass
  generic_single_entry("Residues to color? Single letter code (e.g. DE or de will color Asp/Glu)","A","Color entered residue types!",color_from_string)

#Mutate active chain to entered sequence
default_seq="MAAAA"
def mutate_by_resnum():
  def enter_seq(seq):
    global default_seq
    seq=str(seq).upper()
    seq.replace(" ", "")
    seq_dic={}
    len_seq=len(seq)
    n=0
    nmax=len_seq+1
    aa_dic={'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS','E':'GLU',
    'Q':'GLN','G':'GLY','H':'HIS','I':'ILE','L':'LEU','K':'LYS','M':'MET',
    'F':'PHE','P':'PRO','S':'SER','T':'THR','W':'TRP','Y':'TYR','V':'VAL'}
    clean_seq=''
    while (n<len_seq): 
      if ((seq[n].isalpha() and (seq[n] in aa_dic))):
        clean_seq=clean_seq+seq[n]
      n=n+1
    seq=clean_seq
    default_seq=seq
    len_seq=len(seq)
    n=0
    while (n<len_seq):
      value=aa_dic[seq[n]]
      seq_dic[n+1]=value
      n=n+1
    mol_id=active_residue()[0]
    ch_id=active_residue()[1]
    sn=0
    last_sn=chain_n_residues(ch_id,mol_id)-1
    turn_off_backup(mol_id)
    while (sn<=last_sn):
      res=resname_from_serial_number(mol_id,ch_id,sn)
      seqnum=seqnum_from_serial_number(mol_id,ch_id,sn)
      ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
      if ((res!=seq_dic.get(seqnum)) and ((res in list(aa_dic.values())) or (res=="MSE"))):
        if seq_dic.get(seqnum):
          mutate(mol_id,ch_id,seqnum,ins_id,seq_dic.get(seqnum))
          if seq_dic.get(seqnum)!="PRO":
            delete_residue_sidechain(mol_id,ch_id,seqnum,ins_id,0)
      sn=sn+1
    turn_on_backup(mol_id)
  generic_single_entry("Enter raw amino acid sequence (must be complete!)",
  default_seq,"Mutate active chain to match sequence using PDB numbering", enter_seq)

#Make an alkyl chain of entered length and autofit to map if available
def make_alkyl_chain():
  def make_alkyl_chain_length_n(n):
    smiles_string=int(int(n)+1)*"c"
    new_molecule_by_smiles_string("",smiles_string,force_libcheck=True)
    delete_hydrogens(molecule_number_list()[-1])
    mol_id=molecule_number_list()[-1]
    ch_id="A"
    res_no=1
    ins_code=""
    altloc=""
    new_residue_name="CXC"
    set_residue_name(mol_id,ch_id,res_no,ins_code,new_residue_name)
    delete_atom(mol_id,ch_id,res_no,ins_code," C  ",altloc)
    prodrg_ify(mol_id,ch_id,res_no,ins_code)
    close_molecule(molecule_number_list()[-1])
    close_molecule(molecule_number_list()[-1])
    if imol_refinement_map()!=-1:
      fit_chain_to_map_by_random_jiggle(mol_id,ch_id,1000,0.1)
  generic_single_entry("How many carbons do you want in the chain?",
    "10","Make alkyl chain",make_alkyl_chain_length_n)

#Make new helix (don't fit)
def place_new_helix():
  def place_new_helix_entry(n):
    get_monomer_no_H("ALA")
    mol_id=model_molecule_list()[-1]
    ch_id=chain_ids(mol_id)[0]
    res_no=1
    ins_code=""
    altloc=""
    delete_atom(mol_id,ch_id,1,ins_code," OXT",altloc)
    res_type="auto"
    for i in range(1,int(n)):
      add_terminal_residue_using_phi_psi(mol_id,ch_id,res_no,
      res_type,-57.82,-47)
      res_no=res_no+1
  generic_single_entry("How many residues for helix?",
  "10","Place helix",place_new_helix_entry)
  
#Make new strand (don't fit)
def place_new_strand():
  def place_new_strand_entry(n):
    get_monomer_no_H("ALA")
    mol_id=model_molecule_list()[-1]
    ch_id=chain_ids(mol_id)[0]
    res_no=1
    ins_code=""
    altloc=""
    delete_atom(mol_id,ch_id,1,ins_code," OXT",altloc)
    res_type="auto"
    for i in range(1,int(n)):
      add_terminal_residue_using_phi_psi(mol_id,ch_id,res_no,
      res_type,-139,135)
      res_no=res_no+1
  generic_single_entry("How many residues for strand?",
  "10","Place strand",place_new_strand_entry)

#Make new 3-10 helix (don't fit)
def place_new_3_10_helix():
  def place_new_3_10_helix_entry(n):
    get_monomer_no_H("ALA")
    mol_id=model_molecule_list()[-1]
    ch_id=chain_ids(mol_id)[0]
    res_no=1
    ins_code=""
    altloc=""
    delete_atom(mol_id,ch_id,1,ins_code," OXT",altloc)
    res_type="auto"
    for i in range(1,int(n)):
      add_terminal_residue_using_phi_psi(mol_id,ch_id,res_no,
      res_type,-49,-26)
      res_no=res_no+1
  generic_single_entry("How many residues for 3-10 helix?",
  "10","Place 3-10 helix",place_new_3_10_helix_entry)

#Renumber active segment by active residue
def renumber_seg_by_active_res():
  def renum_seg(new_num):
    new_num=int(new_num)
    current_num=active_residue()[2]
    mol_id=active_residue()[0]
    ch_id=active_residue()[1]
    seg_count=0
    segments=segment_list(mol_id)
    last_res=last_polymer_residue(mol_id,ch_id)
    print(last_res)
    first_res=first_residue(mol_id,ch_id)
    print(first_res)
    new_seg_list=[]
    for seg in segments:
      if ch_id==seg[1]:
        new_seg_list.append(seg) 
    print(new_seg_list)
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
        delete_all_extra_restraints(mol_id)
        set_show_extra_restraints(mol_id,0)
        set_show_extra_restraints(mol_id,1)
  generic_single_entry("New number for this residue?",
  str(active_residue()[2]),"Renumber",renum_seg)

def find_sequence_with_entry():
  generic_single_entry("Enter sequence fragment to find",
  "MAAAA","Find sequence in active chain",find_sequence_in_current_chain)

def user_defined_add_arbitrary_length_bond_restraint(bond_length=2.0):
  def make_restr(text_list, continue_qm):
    s = "Now click on 2 atoms to define the additional bond restraint"
    add_status_bar_text(s)
    dist = text_list[0]
    try:
      bl = float(dist)
    except:
      bl = False
      add_status_bar_text("Must define a number for the bond length")
    if bl:
      def make_restr_dist(*args):
        atom_spec_1 = args[0]
        atom_spec_2 = args[1]
        imol = atom_spec_1[1]
        print("BL DEBUG:: imol: %s spec 1: %s and 2: %s" %(imol, atom_spec_1, atom_spec_2))
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
'''
#
# The full custom menu block is preserved verbatim below. Remove the opening
# and closing triple quotes around the archive to re-enable it.
GTK_UI_ARCHIVE = r'''
#****Make and arrange menus****

if GUI_PYTHON_AVAILABLE:
  menu = attach_module_menu_button("Custom")

  # Coot 1.1 builds can expose coot_gui without a directly importable Gio API.
  # Flatten the menu into one module button for maximum startup compatibility.
  submenu_display = menu
  submenu_fit = menu
  submenu_renumber = menu
  submenu_settings = menu
  submenu_build = menu
  submenu_mutate = menu
  submenu_copy = menu
  submenu_delete = menu
  submenu_merge = menu
  submenu_maps = menu
else:
  menu = None
  submenu_display = submenu_fit = submenu_renumber = submenu_settings = None
  submenu_build = submenu_mutate = submenu_copy = submenu_delete = None
  submenu_merge = submenu_maps = None

#**** Populate submenus ****
#"Display..."

add_simple_coot_menu_menuitem(submenu_display, "All Molecules use \"C-alpha\" Symmetry", lambda func: [valid_model_molecule_qm(imol) and symmetry_as_calphas(imol, 1) for imol in molecule_number_list()])

add_simple_coot_menu_menuitem(submenu_display, "Toggle Symmetry", 
lambda func: set_show_symmetry_master(not get_show_symmetry()))

add_simple_coot_menu_menuitem(submenu_display, "Clear labels and distances", 
lambda func: clear_distances_and_labels())

add_simple_coot_menu_menuitem(submenu_display,
"Switch all mols to CA representation",lambda func: all_mols_to_ca())

add_simple_coot_menu_menuitem(submenu_display,
"Color active mol by rotamer prob (outliers magenta) and missing atoms (blue)", lambda func: color_rotamer_outliers_and_missing_atoms(active_residue()[0]))

add_simple_coot_menu_menuitem(submenu_display,
"Color active mol by hydrophobics (orange), polars (blue), glys (magenta) and pros (green)", lambda func: color_polars_and_hphobs(active_residue()[0]))

add_simple_coot_menu_menuitem(submenu_display,
"Color active mol by charge (+ve blue, -ve red)", lambda func: color_by_charge(active_residue()[0]))

add_simple_coot_menu_menuitem(submenu_display,
"Uncolor other chains in active mol", lambda func: uncolor_other_chains())

add_simple_coot_menu_menuitem(submenu_display,
"Color active chain", lambda func: color_active_chain())

add_simple_coot_menu_menuitem(submenu_display,
"Color active segment", lambda func: colour_active_segment())

add_simple_coot_menu_menuitem(submenu_display,
"Color by protein/nucleic acid", lambda func: color_protein_na(active_residue()[0]))


add_simple_coot_menu_menuitem(submenu_display,
"Color waters", lambda func: color_waters(active_residue()[0]))

add_simple_coot_menu_menuitem(submenu_display, "Colour entered subset of protein residues for active mol", lambda func: color_protein_residue_subset())

add_simple_coot_menu_menuitem(submenu_display,
"Color active mol by CaBLAM outliers (blue) (needs phenix)", lambda func: color_by_cablam2(active_residue()[0]))

add_simple_coot_menu_menuitem(submenu_display,
"Color active mol by EMringer outliers (red) (needs phenix)", lambda func: color_emringer_outliers(active_residue()[0],scroll_wheel_map()))

add_simple_coot_menu_menuitem(submenu_display,
"Color active mol by ramachandran outliers (blue) (needs phenix)", lambda func: color_by_rama(active_residue()[0]))

add_simple_coot_menu_menuitem(submenu_display,
"Color active mol by map/model-CC at 4 A (Slow, needs phenix, P1 only)", lambda func: color_by_cc(active_residue()[0]))

add_simple_coot_menu_menuitem(submenu_display, "Highlight chain breaks in active mol", lambda func: highlight_chain_breaks())

add_simple_coot_menu_menuitem(submenu_display, "Highlight chain breaks in all mols", lambda func: highlight_all_chain_breaks())

add_simple_coot_menu_menuitem(submenu_display, "Open current view in UCSF Chimera. Experimental!", lambda func: open_in_chimera())

add_simple_coot_menu_menuitem(submenu_display, "Show local probe dots (H-bonds, VdW and baddies only)", lambda func: local_hbonds_and_baddies())

add_simple_coot_menu_menuitem(submenu_display, "Show global probe dots (bad overlaps only)", lambda func: global_hbonds_and_baddies())

add_simple_coot_menu_menuitem(submenu_display, "Clear probe dots", lambda func: clear_dots())

add_simple_coot_menu_menuitem(submenu_display, "Find sequence in active chain", lambda func: find_sequence_with_entry())

add_simple_coot_menu_menuitem(submenu_display, "Make stereo-wiggle GIF of current scenre (Needs ImageMagick!)", lambda func: make_rot_gif())



#"Fit..."
add_simple_coot_menu_menuitem(submenu_fit, "Fit all chains to map", 
lambda func: rigid_fit_all_chains())

add_simple_coot_menu_menuitem(submenu_fit, "Stepped sphere refine active chain",
lambda func: stepped_sphere_refine(active_residue()[0],active_residue()[1]))

add_simple_coot_menu_menuitem(submenu_fit, "Fit current chain to map", 
lambda func: rigid_fit_active_chain())

add_simple_coot_menu_menuitem(submenu_fit, 
"Jiggle-fit current chain to map (Slow!)", lambda func: jiggle_fit_active_chain())

add_simple_coot_menu_menuitem(submenu_fit,
"Jiggle-fit current chain to B-smoothed map (Slow!)", lambda func: jiggle_fit_active_chain_smooth())

add_simple_coot_menu_menuitem(submenu_fit, "Jiggle-fit all chains to map (very slow!)",
lambda func: jiggle_fit_all_chains())

add_simple_coot_menu_menuitem(submenu_fit, 
"Jiggle-fit current mol to map (Slow!)", lambda func: jiggle_fit_active_mol())

add_simple_coot_menu_menuitem(submenu_fit, 
"Fit polyala loop (click start and end)", lambda func: fit_polyala_gui())

add_simple_coot_menu_menuitem(submenu_fit, "Fit all segments", lambda func: rigid_body_fit_segments())

add_simple_coot_menu_menuitem(submenu_fit, "Fit this segment", lambda func: fit_this_segment())

add_simple_coot_menu_menuitem(submenu_fit,
"Prosmart self restrain active mol",lambda func: run_prosmart_self())

add_simple_coot_menu_menuitem(submenu_fit,
"Cylinder refine (click start and end of range)",lambda func: refine_residues_sphere_click())

add_simple_coot_menu_menuitem(submenu_fit, "Add distance restraint (click two atoms)",lambda func:  user_defined_add_arbitrary_length_bond_restraint(bond_length=2.0))

#"Renumber..."

add_simple_coot_menu_menuitem(submenu_renumber, "Renumber active chain by first res", 
lambda func: renumber_by_first_res())

add_simple_coot_menu_menuitem(submenu_renumber, 
"Renumber active chain by last res", lambda func: renumber_by_last_res())

add_simple_coot_menu_menuitem(submenu_renumber,
"Renumber active chain by current res", lambda func: renumber_by_active_res())

add_simple_coot_menu_menuitem(submenu_renumber,"Renumber from N-term to active residue",
lambda func: renumber_n_term_segment())

add_simple_coot_menu_menuitem(submenu_renumber,"Renumber from active residue to C-term",
lambda func: renumber_c_term_segment())

add_simple_coot_menu_menuitem(submenu_renumber, "Renumber segment by active res", lambda func: renumber_seg_by_active_res())


#"Settings..."

add_simple_coot_menu_menuitem(submenu_settings,
"Auto-scale B-factor coloring for active mol",lambda func: autoscale_b_factor())

add_simple_coot_menu_menuitem(submenu_settings,
"Set Bfac for new atoms to mean B for active mol",lambda func: set_new_atom_b_fac_to_mean())


#"Build..."
add_simple_coot_menu_menuitem(submenu_build,
"Forced addition of terminal residue (click terminus)", 
lambda func: force_add_terminal_residue())

add_simple_coot_menu_menuitem(submenu_build,
"Grow helix (click terminus)",
lambda func: grow_helix())

add_simple_coot_menu_menuitem(submenu_build,
"Grow strand (click terminus)",
lambda func: grow_strand())

add_simple_coot_menu_menuitem(submenu_build,
"Grow parallel strand (click terminus)",
lambda func: grow_parallel_strand())

add_simple_coot_menu_menuitem(submenu_build,
"Grow 3-10 helix (click terminus)",
lambda func: grow_helix_3_10())

add_simple_coot_menu_menuitem(submenu_build, 
"Shorten loop by one residue", lambda func: shorten_loop())

add_simple_coot_menu_menuitem(submenu_build, 
"Lengthen loop by one residue", lambda func: lengthen_loop())

add_simple_coot_menu_menuitem(submenu_build,
"Get fractional coordinates of active atom",lambda func: get_fract_coords()) 

add_simple_coot_menu_menuitem(submenu_build,
"Common monomers",lambda func: pick_common_monomers())

add_simple_coot_menu_menuitem(submenu_build, "Make alkyl chain of length n", lambda func: make_alkyl_chain())

add_simple_coot_menu_menuitem(submenu_build, "Make alpha helix of length n", lambda func: place_new_helix()) 

add_simple_coot_menu_menuitem(submenu_build, "Make 3-10 helix of length n", lambda func: place_new_3_10_helix())

add_simple_coot_menu_menuitem(submenu_build, "Rebuild backbone (click start,end)", lambda func: rebuild_backbone_wrapper())

add_simple_coot_menu_menuitem(submenu_build, "Rebuild and reverse backbone of segment (click segment)", lambda func: rebuild_backbone_reverse_wrapper())





#"Mutate...

add_simple_coot_menu_menuitem(submenu_mutate,
"Mutate range to UNK (click start and end)", lambda func: mutate_residue_range_by_click_a())

add_simple_coot_menu_menuitem(submenu_mutate,
"Mutate range to ALA (click start and end)", lambda func: mutate_residue_range_by_click_ala_a())

add_simple_coot_menu_menuitem(submenu_mutate, "Mutate all Mets to MSE", lambda func: mutate_all_mets_to_mse())

add_simple_coot_menu_menuitem(submenu_mutate, "Mutate all MSEs to Met", lambda func: mutate_all_mse_to_met())

add_simple_coot_menu_menuitem(submenu_mutate,
"Mutate active chain to template sequence (numbering must match sequence!)", lambda func: mutate_by_resnum())

#"Copy..."
add_simple_coot_menu_menuitem(submenu_copy, "Copy current chain", 
lambda func: copy_active_chain())

add_simple_coot_menu_menuitem(submenu_copy, "Cut current chain", 
lambda func: cut_active_chain())

add_simple_coot_menu_menuitem(submenu_copy, "Copy active segment", lambda func: copy_active_segment())

add_simple_coot_menu_menuitem(submenu_copy, "Cut active segment", lambda func: cut_active_segment())

add_simple_coot_menu_menuitem(submenu_copy,
"Copy fragment (click start and end)", lambda func: copy_frag_by_click())

add_simple_coot_menu_menuitem(submenu_copy,
"Cut fragment (click start and end)", lambda func: cut_frag_by_click())

add_simple_coot_menu_menuitem(submenu_copy,
"Copy active chain to NCS equivs", lambda func: copy_ncs_chain_from_active())


#"Delete..."
add_simple_coot_menu_menuitem(submenu_delete,
"Delete active chain", lambda func: delete_chain())

add_simple_coot_menu_menuitem(submenu_delete, "Delete active segment", lambda func: delete_active_segment())


add_simple_coot_menu_menuitem(submenu_delete, 
"Delete hydrogens from molecule", lambda func: delete_h_active())

add_simple_coot_menu_menuitem(submenu_delete,
"Delete sidechains in range (click start and end)", lambda func: delete_sidechain_range_by_click_a())  

#"Merge..."

add_simple_coot_menu_menuitem(submenu_merge, 
"Merge two mols (click two; 2nd into 1st)", lambda func: merge_fragments())

add_simple_coot_menu_menuitem(submenu_merge, "Merge chains (click two; 2nd into 1st)", lambda func: merge_chains())


#"Maps..."
add_simple_coot_menu_menuitem(submenu_maps,
"Sharpen (enter B-factor)",lambda func: sharpen_by_entered_factor())

add_simple_coot_menu_menuitem(submenu_maps,
"Change hi-res limit for map",lambda func: change_hires_limit_copy())

add_simple_coot_menu_menuitem(submenu_maps,
"Go to center of scrollable map",lambda func: goto_center_of_map())

add_simple_coot_menu_menuitem(submenu_maps,
"Set refinement map to scrollable map",lambda func: set_map_to_scrollable_map())
'''
