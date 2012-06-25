from fitcmp import VERSION
import tarfile, os, glob, sys

scriptpath = os.path.abspath(os.path.dirname(sys.argv[0]))

archive_path = "%s/fitcmp_%.2f.tar.gz" % (scriptpath, VERSION)
files = []
files.append("%s/fitcmp.py" % scriptpath)
files.extend(glob.glob("%s/*.txt" % scriptpath))
files.extend(glob.glob("%s/*.cfg" % scriptpath))
files.extend(glob.glob("%s/vpstuff/*.py" % scriptpath))

with tarfile.open(archive_path, "w:gz") as tar:
    for f in files:
        tar.add(f, arcname=f[f.index('fitcmp'):])
