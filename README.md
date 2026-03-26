# coot1-trimmings - Customizations & scripts for Coot 1.x

`coot1_trimmings.py` is a Coot startup script for Coot 1.1.x, ported from `coot-trimmings` for Coot 0.9.x with the assistence of the Codex LLM. It adds custom keybindings, map-display presets, model-editing shortcuts, and a small set of startup preferences.

Custom items that rely on Gtk menus and dialogs are not working yet - stay tuned! **EDIT:** Just added them - docs to come!

A quick demo of some of the keyboard actions is here:

https://www.dropbox.com/scl/fi/1wytsp29147x5yt83coz9/coot1_trimmings_demo.mov?rlkey=7gjcty91ucuphhcoml14863gy&dl=0

And a demonstration of keyboard copy/paste of ligand/solvent molecules is here:

https://www.dropbox.com/scl/fi/fj4nod1axb4vhvz0tewtc/smart_copy_paste.mov?rlkey=b8p7funnrexxfojytlww1obde&dl=0

_Note: if you are looking for trimmings for Coot 0.9, you are in the wrong place - please go here instead: https://github.com/olibclarke/coot-trimmings_

## Install

Copy ￼`coot1_trimmings.py`￼ into Coot's startup-script directory (e.g. ~/.config/Coot/ on Mac) and restart Coot.

Example:

```bash
cp coot1_trimmings.py ~/.config/Coot/
```
(I don't know the location on Linux or Windows, but it will be wherever your `coot_preferences.py` file for Coot 1 resides.)

_Note: If you see an errror like <!DOCTYPE html> ^ SyntaxError: invalid syntax, you have downloaded the webpage rather than the actual Python script. Make sure you clicked "download raw file" when downloading._

## Startup behavior

The script changes a few defaults at startup, including:

- Left mouse button used for view rotation
- Brighter symmetry molecule colours
- Disabled dragged map and smooth recentering.

## Example use case/Tutorial
1. Load up an EM-map and model (e.g. 8HEZ & EMD-34705)
2. Adjust the contour - `"Shift-1"` - `"Shift-9"` set the level from 1 to 9 sigma respectively, and `"+"` and `"_"` nudge the contour level up and down. Try adjusting the map radius - `;` to decrease, `'` to increase.
3. Zoom in - you might notice the mesh is a bit coarse, which can make fine details difficult to interpret. This is a 2.8Å structure, but was solved at 1.1Å/pixel. Let's try resampling and brightening things up a bit to make it easier to see. to do this press the `"C"` shortcut key. (_Note: there are some smarts here - if the map is an X-ray map, it will not be resampled, but color changes etc will still happen; if it is an EM map but is already finer than 0.5 Å/pix, the sampling will not be changed; If it is coarser than 0.5Å/pix, it will be resampled to 0.5Å/pix and the original map will be closed.)_
4. Try adjusting the clipping - use `"-"` to increase the clipping (making the slab tighter around the center of rotation), and `"="` to make it wider. Use `";"` and `" ' "` to adjust the radius of the map in local mode.
5. The local mesh is great for looking at a local region - say a few sidechains - but what if we want to see features of the global map? To switch to a surface view of the whole map, press `"G"`. Press `"G"` again to switch back to the local view. When in global mode, you can also use `":"` and `" " "` to adjust the opacity. (_Note: On my Mac Coot 1 stil has some graphical glitches when adjusting opacity, I suspect because of Gtk4 issues_)
6. Sometimes, you want to quickly toggle the model off to better see the map, or vice-versa. To toggle the model on/off, press `"/"`; to switch between multiple models, press `"Shift-/"`. To toggle the map, press `"\`"`; To switch between multiple maps, press `"Shift-\`"`.
7. Try cycling the representation of the model - press the `"\["` and `"\]"` keys to move back and forward between the different modes (e.g. CA-only and all atom).
8. Try cylinder refine - center on an atom, and press `"a"`. It will refine +/- 5 residues around the center, and all segments within 4Å of the primary range.
9. Try quick refine zone - press `"A"`, then click two atoms to define a range and refine.
10. Try copy/pasting a ligand - center on a ligand or solvent molecule, and press `"c"`. A status message in the lower bar should indicate that the molecule has been copied. Move somewhere else, and press `"v"`. A copy should be pasted and merged into the active molecule. You can repeat this operation multiple times.

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
| `+` | Increase current map contour by 0.1 sigma |
| `_` | Decrease current map contour by 0.1 sigma |
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
- _Caveat: on Mac there is currently a bug where after applying self restraints (either custom or native), applying Keyboard Mutate will cause Coot to crash. This is the case on 1.1.20, might be fixed in daily build._

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

