-- Root module for this repository's standalone Lean build.
-- It is NOT a copy: the shared project's root also imports other, unrelated targets.
-- Here it imports only the Hart–Rowley files, which themselves import no library.
-- `lean-toolchain`, `lakefile.toml` and `lake-manifest.json` are unmodified copies
-- from the shared project (Lean v4.37.0-rc1; Mathlib pinned in the manifest but not imported).
import Pipeline.HartRowley.Main
