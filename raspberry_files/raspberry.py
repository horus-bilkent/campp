import subprocess
import time
import random
import jsonpickle
import json
import os.path
from pykafka import KafkaClient
  
CONFIG_FILE = './raspberry_config.json'

def generate_initialization_message(appliance_id, reply_to, image_1, image_2):
	with open(image_1, 'rb') as image_file:
		image_1_bytes = bytearray(image_file.read())
	
	with open(image_2, 'rb') as image_file:
		image_2_bytes = bytearray(image_file.read())
		
	message = {}
	message['appliance_id'] = appliance_id
	message['image_1'] = jsonpickle.encode(image_1_bytes)
	message['image_2'] = jsonpickle.encode(image_2_bytes)
	message['reply_to'] = reply_to
	message['timestamp'] = str(int(time.time()))
	return json.dumps(message)

def generate_stream_message(appliance_id, value):
	message = {}
	message['appliance_id'] = appliance_id
	message['value'] = str(value)
	message['timestamp'] = str(int(time.time()))
	return json.dumps(message)

if __name__ == "__main__":
	with open(CONFIG_FILE) as config_file:    
		config = json.load(config_file)
		
	try:
		
		broker_connect = config['broker_connect']
		appliance_id = config['appliance_id']
		initialization_topic = config['initialization_topic']
		read_topic = config['read_topic']
		write_initialization = config['write_initialization']
		find_needle = config['find_needle']
		diff_images = config['diff_images']
		granularity = config['read_granularity']
		stream_topic = config['stream_topic']
	except KeyError, e:
		print 'Error in config file: ' + str(e)
		
	client = KafkaClient(hosts=broker_connect)
	topic_initialization = client.topics[bytes(initialization_topic)]
	topic_read = client.topics[bytes(read_topic)]
	topic_stream = client.topics[bytes(stream_topic)]
	if not os.path.isfile(write_initialization):
		# image_1 = './image_1_' + str(int(time.time()))
		# image_2 = './image_2_' + str(int(time.time()))
		# ret_val = subprocess.call([diff_images, image_1, image_2])
		
		# to try
		image_1 = './image_1.jpg'
		image_2 = './image_2.jpg'
		ret_val = 0
		
		if  ret_val != 0:
			print 'Initialization, Find two different images failed.'
			exit(1)
		
		else:
			with topic_initialization.get_sync_producer(max_request_size=50000000) as producer:
				# generate the json message
				print 'Sending'
				message = generate_initialization_message(appliance_id, read_topic, image_1, image_2)
				producer.produce(message)
		
		# read the initialization parameters
		consumer = topic_read.get_simple_consumer()
		for message_recv in consumer:
			if message_recv is not None:
				message_recv_json = json.loads(message_recv.value)
				if 'timestamp' in message_recv_json and int(time.time()) <= int(message_recv_json['timestamp']) + 15:
					print 'Received'
					with open(write_initialization, 'w') as outfile:
						outfile.write(message_recv_json['value'])
						break
				else:
					print 'Skipping old message'
	else:
		print 'Initialization config file already exists, skipping initialization.'

	
	# after initialization call findNeedle
	with topic_stream.get_sync_producer() as producer:
			child = subprocess.Popen([find_needle, write_initialization], shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE)		
			while True:
				child.stdin.write("Send\n")
				child.stdin.flush()
				line = child.stdout.readline()
				child.stdout.flush()
				value = float(line)
				message = generate_stream_message(appliance_id, value)
				print 'Sending...'
				print message
				producer.produce(message)
				time.sleep(granularity)
				
