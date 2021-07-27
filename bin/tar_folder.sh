#!/bin/bash

destination=$1
folder=$2

tar -zcf $destination $folder

if [ ! $? -e 0 ]
then
	echo "Command returned non-zero return value. File should be checked. $destination"
fi

