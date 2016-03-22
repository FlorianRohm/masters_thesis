#!/bin/zsh

SCRIPT=$(readlink -f $0)
SCRIPTPATH=`dirname $SCRIPT`

docker run -it -p 8888:8888  -v $SCRIPTPATH:/opt/notebooks/mine walberla/runenv-ubuntu-python
