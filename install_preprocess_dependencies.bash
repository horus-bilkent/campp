sudo -v
if [ $? -neq 0 ]
then
	echo 'RUN WITH SUDO!'
	exit 255
fi

git clone -b master https://github.com/cedricve/raspicam

cd raspicam

mkdir build

cd build

cmake ..

make

make install

ldconfig

rm -rf ../../raspicam

pip install jsonpickle

pip install pykafka