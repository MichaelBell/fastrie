xptMiner
========

This is a Riecoin (RIC)-focused release for the Parallella board,
built upon xptMiner.  It contains only a Riecoin miner - for the 
reference implementation of xptMiner, please see 
https://github.com/jh000/xptMiner

Authors:  
 * xptMiner was written by jh00;
 * This version is based upon dga's fastrie variant, which in turn
   was based upon the Unix port by Clintar
 * The Parallella implementation is by Michael Bell, aka rockhawk

Some instructions to get started:

PREREQUISITES 
=============
Linaro:

    sudo apt-get install build-essential m4 openssl libssl-dev git libjson0 libjson0-dev libcurl4-openssl-dev wget libgmp-dev

POOL SET UP
===========

Create an account on http://ypool.net and configure a name and password for
a Riecoin worker.  

BUILDING
========

    git clone https://github.com/MichaelBell/fastrie.git
    cd fastrie/xptMiner
    make -j2

RUNNING
=======

    ./run.sh -u username.riecoinworkername -p workerpassword

This has a default 2% donation that can be set using the -d option (-d 2.5 would be 2.5% donation)

RUNNING AS A SERVICE
====================

To dedicate a Parallella to running the miner, edit rieminer.conf in the
fastrie directory to specify your own username and password, then
copy it to /etc/init.  After that

    sudo service rieminer start

will begin mining.

Note this service will automatically reboot the Parallella to recover on error,
after a one minute sleep so you can stop the service if it is repeatedly 
rebooting.  Stop it with:

    sudo service rieminer stop

