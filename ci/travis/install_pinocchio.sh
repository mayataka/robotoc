sudo apt install -qqy lsb-release gnupg2 curl
echo "deb [arch=amd64] http://robotpkg.openrobots.org/packages/debian/pub $(lsb_release -cs) robotpkg" | sudo tee /etc/apt/sources.list.d/robotpkg.list
curl http://robotpkg.openrobots.org/packages/debian/robotpkg.key | sudo apt-key add 
sudo apt update
sudo apt install -qqy robotpkg-py27-pinocchio
echo export PATH=/opt/openrobots/bin:$PATH
echo export PKG_CONFIG_PATH=/opt/openrobots/lib/pkgconfig:$PKG_CONFIG_PATH
echo export LD_LIBRARY_PATH=/opt/openrobots/lib:$LD_LIBRARY_PATH
echo export PYTHONPATH=/opt/openrobots/lib/python2.7/site-packages:$PYTHONPATH
echo export CMAKE_PREFIX_PATH=/opt/openrobots:$CMAKE_PREFIX_PATH