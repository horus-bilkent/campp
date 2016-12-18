import subprocess
import time
import random
import jsonpickle
import json
import os.path
import os
import sys
import picamera
import auxiliary
from pykafka import KafkaClient

def capture_image(image_path, image_res):
		camera = PiCamera()
		camera.resolution = image_res
		sleep(2)
		camera.capture(image_path)


CONFIG_FILE = './raspberry_config.json'

if __name__ == "__main__":
	if sys.argc == 1:
		config = auxiliary.read_config(CONFIG_FILE)
	else:
		config = auxiliary.read_config(sys.argv[1])
			
	try:
		
		broker_connect = config['broker_connect']
		appliance_id = config['appliance_id']
		appliance_topic = config['appliance_topic']
		write_initialization = config['write_initialization']
		resolution = eval(config['capture_resolution'])
		find_needle = config['find_needle']
		granularity = config['read_granularity']
		stream_topic = config['stream_topic']
	except KeyError, e:
		print 'Error in config file: ' + str(e)
		
	client = KafkaClient(hosts=broker_connect)
	topic_appliance = client.topics[bytes(read_topic)]
	
	if not os.path.isfile(write_initialization):
		consumer = topic_appliance.get_simple_consumer()
		for message_recv in consumer:
			print 'Waiting for client initialization...'
			if message_recv is not None:
				begin_msg = json.loads(message_recv.value)
				if 'transaction_id' in begin_msg and 'type' in begin_msg and begin_msg['type'] == 'initialization_begin':
					print 'Received, initialization command.'
					transaction_id = begin_msg['transaction_id']

		# capture image
		image_path = './image_ ' + auxiliary.get_time() + '.jpg'
		capture_image(image_path, resolution)
		with open(image_path, 'rb') as image_file:
			image_pickled = jsonpickle.encode(bytearray(image_file.read()))
		
		with topic_initialization.get_sync_producer(max_request_size=50000000) as producer:
			# generate the json message
			print 'Sending'
			image_message = auxiliary.generate_message(appliance_id=appliance_id, msg_type='initialization_image', transaction_id=transaction_id, value=image_pickled)
			producer.produce(image_message)
			# os.remove(image)

		# read the initialization parameters after client initialization
		consumer = topic_read.get_simple_consumer()
		for message_recv in consumer:
			if message_recv is not None:
				end_msg = json.loads(message_recv.value)
				if 'transaction_id' in end_msg and 'type' in end_msg and end_msg['type'] == 'initialization_end' and end_msg['transaction_id'] == transaction_id:
					print 'Received'
					with open(write_initialization, 'w') as outfile:
						outfile.write(message_recv_json['value'])
						break
		
	else:
		print 'Initialization file already exists, skipping initialization.'
	
	# after initialization call findNeedle continously
	with topic_stream.get_sync_producer() as producer:
			child = subprocess.Popen([find_needle, write_initialization], shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE)		
			while True:
				image_path = 'image_' + auxiliary.get_time() + '.jpg'
				capture_image(image, resolution)
				print 'Processing: ' + image_path
				child.stdin.write(image_path)
				child.stdin.flush()
				line = child.stdout.readline()
				child.stdout.flush()
				value = float(line)
				message = auxiliary.generate_message(appliance_id=appliance_id, value=value)
				print 'Sending...'
				print message
				producer.produce(message)
				os.remove(image_path)
				time.sleep(granularity)
