import os
import subprocess
from .MK_read.main import dump_dat
from pathlib import Path

def convert(tar_file,typelabel_dat):
    """
    Get a typelabel-only (for usage with a KIM SM) data file from a CHARMM-GUI .tgz file
    Args:
        tar_file: CHARMM-GUI tgz
        typelabel_dat: path to data file to write. Parent dirs will be created if they don't exist
    """
    subprocess.check_output("tar --wildcards --strip-components=2 -xkf %s charmm-gui*/lammps/step3_input.data"%tar_file,shell=True)    
    Path(os.path.dirname(typelabel_dat)).mkdir(parents=True, exist_ok=True)    
    dump_dat('step3_input.data',typelabel_dat)
    os.remove('step3_input.data')

