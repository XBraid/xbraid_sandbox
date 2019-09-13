# To Run, 
# export PYTHONPATH="/home/jbschroder/joint_repos/braid_sandbox/xbraid-with-cython/cyclobraid/examples/Cyclops-Lite-master/source:$HOME/.local/lib/python3.6"
# ipython 3
# >>> import cyclobraid
# >>> cyclobraid.braid_init_py()

rm -rf build
python3 setup-cyclobraid.py install --prefix=$HOME/.local

