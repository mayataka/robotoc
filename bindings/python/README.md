## Configure environment variables for Python binding
The Python binding will be installed in `IDOCP_INSTALL_DIR/lib/python3.8/site-packages` where `IDOCP_INSTALL_DIR` is the install directory of `idocp` configured in CMake (e.g., by `-DCMAKE_INSTALL_PREFIX`).
To use the installed Python library, it is convenient to set the environment variable as

```
export PYTHONPATH=IDOCP_INSTALL_DIR/lib/python3.8/site-packages:$PYTHONPATH 
```

Note that if you use another Python version than `python3.8`, please adapt it.