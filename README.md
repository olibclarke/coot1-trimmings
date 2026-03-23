# coot1-trimmings - Customizations & scripts for Coot 1.x

`coot1_trimmings.py` is a Coot startup script for Coot 1.1.x, ported from `coot-trimmings` for Coot 0.9.x with the assistence of the Codex LLM. It adds custom keybindings, map-display presets, model-editing shortcuts, and a small set of startup preferences.

Custom items that rely on Gtk menus and dialogs are not working yet - stay tuned!

## Install

Copy ￼`coot1_trimmings.py`￼ into Coot's startup-script directory (e.g. ~/.config/Coot/ on Mac) and restart Coot.

Example:

```bash
cp coot1_trimmings.py ~/.config/Coot/
```

## Startup behavior

The script changes a few defaults at startup, including:

- Left mouse button used for view rotation
- Mouse wheel map scrolling disabled (doesn't work well on Mac trackpad) - you can re-enable by changing `set_scroll_by_wheel_mouse(0)` to `set_scroll_by_wheel_mouse(1)` in the script. Otherwise, I have included shortcuts described below for keyboard adjustment of contour.
- Brighter symmetry molecule colours
- Disabled dragged map and smooth recentering.

## Custom keybindings

These are the main custom bindings currently defined by the script.

***Notes:***

- Several bindings override native Coot 1 shortcuts (documented below).
- Most commands act on the active residue, active model, or active scrollable map.

### Display and navigation

| Key     | Action                                                       |
|---------|--------------------------------------------------------------|
| `G`     | Toggle the active map between a local mesh view and a global solid-surface view, expanding map radius and clipping to fit the whole map and restoring them on return |
| `?`     | Show only the active model, or cycle displayed models if there are multiple displayed |
| `~`     | Show only the active map, or cycle displayed maps if there are multiple displayed |
| `` ` `` | Toggle display of all maps                                   |
| `/`     | Toggle display of all model molecules                        |
| `[`     | Cycle model representation mode forward                      |
| `]`     | Cycle model representation mode backward                     |
| `{`     | Cycle symmetry representation mode forward                   |
| `}`     | Cycle symmetry representation mode backward                  |
| `>`     | Go to next residue in chain                                  |
| `<`     | Go to previous residue in chain                              |
| `O`     | Go to the equivalent residue on the NCS master chain         |
| `D`     | Toggle environment distances / H-bond display                |
| `Z`     | Clear distances and labels                                   |

### Map controls

| Key | Action |
| --- | --- |
| `!` | Set current map to 1 sigma |
| `@` | Set current map to 2 sigma |
| `#` | Set current map to 3 sigma |
| `$` | Set current map to 4 sigma |
| `%` | Set current map to 5 sigma |
| `^` | Set current map to 6 sigma |
| `&` | Set current map to 7 sigma |
| `*` | Set current map to 8 sigma |
| `(` | Set current map to 9 sigma |
| `|` | Increase current map contour by 0.5 sigma |
| `_` | Decrease current map contour by 0.5 sigma |
| `'` | Increase map radius |
| `;` | Decrease map radius |
| `"` | Increase active-map surface opacity in global solid mode |
| `:` | Decrease active-map surface opacity in global solid mode |
| `-` | Narrow clipping slab symmetrically about the rotation centre |
| `=` | Widen clipping slab symmetrically about the rotation centre |

Notes:

- `"` and `:` only act when the active map is in the script's global solid-surface mode.
- The preferred solid-surface opacity is remembered per map and restored when switching back into global view.

### Building and editing

| Key | Action                                                       |
|-----|--------------------------------------------------------------|
| `h` | Place a helix                                                |
| `m` | Measure distance (Click two atoms)                           |
| `b` | Go to the nearest density peak around the current rotation centre, using the active map's current contour sigma |
| `w` | Place water at pointer                                       |
| `W` | Place water and refine                                       |
| `y` | Add terminal residue                                         |
| `Y` | Cycle terminal-residue phi                                   |
| `T` | Cycle terminal-residue psi                                   |
| `X` | Delete active residue                                        |
| `K` | Delete active sidechain                                      |
| `k` | Fill partial sidechain                                       |
| `q` | Pepflip active residue                                       |
| `c` | Smart-copy the active ligand, ion, or water                  |
| `v` | Smart-paste the last copied ligand, ion, or water at the pointer |
| `C` | Change the default appearance of the active map; if it is an EM map, resample it to 0.5Å/pixel unless it is already sampled more finely. |

Notes:

- Smart copy only copies non-polymer residues.
- Smart paste can be repeated multiple times from the same copied residue.
- Smart copy currently refuses insertion-code residues instead of copying them ambiguously.
- `b` searches locally around the current rotation centre in a radius of 8 Å.

### Refinement and fitting

| Key | Action                                                       |
|-----|--------------------------------------------------------------|
| `g` | Generate smart local extra Geman-McClure distance restraints for the active model |
| `A` | Refine a clicked residue range (Click two atoms)             |
| `a` | Start a local cylinder refinement around the active residue  |
| `r` | Refine three residues centred on the active residue          |
| `J` | Jiggle-fit the active residue (excludes polymer residues)    |
| `R` | Cycle rotamers for the active residue                        |

Notes:

- `g` clears existing extra restraints on the active model before rebuilding them.
- It adds only inter-residue interatomic `GM` restraints for pairs closer than `3.7 A`.
- Restraints between residues with a gap of more than 10 residues are pruned.
- Long-range backbone `N···O` hydrogen-bond restraints are preserved, including inter-strand and inter-chain contacts.
- It excludes same-residue contacts.
- The idea is to maintain local geometry restraints within secondary structural elements, while allowing them to move around with respect to one another.

### Saving and history

| Key | Action                                                       |
|-----|--------------------------------------------------------------|
| `Q` | Save and overwrite the active model in the current directory, with backup |
| `z` | Undo on active model                                         |
| `x` | Redo on active model                                         |


## Native Coot 1 shortcuts

The following is a (possibly incomplete) list of built-in Coot 1 shortcuts (excluding those that are overwritten by `coot1-trimmings`).

### Single-key shortcuts

| Key      | Native Coot action                      |
|----------|-----------------------------------------|
| `d`      | Step right                              |
| `e`      | Auto-fit rotamer                        |
| `f`      | Decrease clipping                       |
| `i`      | Spin                                    |
| `j`      | Auto-fit rotamer                        |
| `o`      | Go to next NCS chain                    |
| `l`      | Label/unlabel active atom               |
| `n`      | Zoom out                                |
| `p`      | Update go-to atom from current position |
| `s`      | Move backward                           |
| `u`      | Undo move                               |
| `M`      | Keyboard mutate                         |
| `Return` | Accept moving atoms                     |
| `Escape` | Reject moving atoms                     |
| `1`      | Expand front clipping plane             |
| `2`      | Reduce front clipping plane             |
| `3`      | Reduce back clipping plane              |
| `4`      | Expand back clipping plane              |

### Ctrl shortcuts

| Key      | Native Coot action                    |
|----------|---------------------------------------|
| `Ctrl-C` | Duplicate model molecule              |
| `Ctrl-D` | Delete residue                        |
| `Ctrl-G` | Show go-to-residue keyboarding window |
| `Ctrl-E` | Eigen-flip residue                    |
| `Ctrl-I` | Residue info                          |
| `Ctrl-L` | Go to ligand                          |
| `Ctrl-W` | Add water                             |
| `Ctrl-Y` | Redo                                  |
| `Ctrl-Z` | Undo                                  |

## Notes

- Many commands assume an active residue.
- Refinement and fitting commands assume an active refinement map.
- Map-navigation and map-display commands assume an active scrollable map.

