#!/bin/zsh

SCRIPT=$(readlink -f $0)
SCRIPTPATH=`dirname $SCRIPT`

docker run -it -v $SCRIPTPATH:/opt/notebooks/mine walberla/runenv-ubuntu-python /bin/bash
