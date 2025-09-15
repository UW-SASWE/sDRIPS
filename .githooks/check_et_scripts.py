import subprocess
import sys

SCRIPTS = ["src/sdrips/evapotranspiration.py", "src/sdrips/sen_corr_evapotranspiration.py", "src/sdrips/stn_corr_evapotranspiration.py"]

# Get staged files
changed = subprocess.check_output(
    ["git", "diff", "--cached", "--name-only"], text=True
).splitlines()

# Checks if any of the scripts were modified
modified = [s for s in SCRIPTS if s in changed]

if modified and len(modified) < len(SCRIPTS):
    print("[WARNING] sDRIPS developer you modified:", ", ".join(modified))
    print(f"Check if the same change is needed in the other scripts: {', '.join(SCRIPTS)}")
    sys.exit(1)  # block the commit until confirmed

sys.exit(0)  # allow commit
