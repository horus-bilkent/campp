#!/usr/bin/env python
import time
import json
import random
import sys
import jsonpickle
import os
import subprocess
import auxiliary

from pykafka import KafkaClient, SslConfig

CONFIG_FILE='./client_config.json'

def write_jpg(jpg_bytes, to_file):
	with open(to_file, 'w') as outfile:
		outfile.write(jsonpickle.decode(jpg_bytes))
	
if __name__ == "__main__":	
	if len(sys.argv) == 1:
		config = auxiliary.read_config(CONFIG_FILE)
	else:
		config = auxiliary.read_config(sys.argv[1])
			
	try:
		
		broker_connect = config['broker_connect']
		preprocess = config['preprocess']
		appliance_topic = config['appliance_topic']
		appliance_id = config['appliance_id']
	except KeyError, e:
		print 'Error in config file: ' + str(e)
		exit(1)
		
	
	print 'Initializing appliance: ' + str(appliance_id) 
	# config = SslConfig(cafile=ssl_ca_cert, certfile=ssl_client_cert, keyfile=ssl_client_private_key_file, password=ssl_key_password)
	# client = KafkaClient(hosts=broker_connect, ssl_config=config)
	client = KafkaClient(hosts=broker_connect)
	
	topic_appliance = client.topics[bytes(appliance_topic)]
	
	transaction_id=random.randint(1, 2^30)
	# inititialize the sequence
	with topic_appliance.get_sync_producer() as producer:
			begin_msg = auxiliary.generate_message(appliance_id=appliance_id, transaction_id=transaction_id, msg_type='initialization_begin')
			producer.produce(begin_msg)
			
	# consume the message with image
	consumer = topic_appliance.get_simple_consumer(fetch_message_max_bytes=50000000)
	for message_recv in consumer:
		if message_recv is not None:
			image_msg = json.loads(message_recv.value)
			if 'transaction_id' in image_msg and 'type' in image_msg and image_msg['type'] == 'initialization_image' and image_msg['transaction_id'] == transaction_id: 
				print 'Received, initialization image.'
				image_bytes = image_msg['value']
				image_path = './image_' + str(int(time.time())) + '.jpg'
				output_file = './preprocess_output_' + str(int(time.time()))
				
				write_jpg(image_bytes, image_path)
				ret_val = 0
				# ret_val = subprocess.call([preprocess, image_name, output_file])
				
				if ret_val != 0:
					print 'Initialization failed, cannot execute preprocessing'
					exit(1)
												
				with open(output_file, 'r') as infile:
					output = infile.read()
					
				os.remove(image_path)
				# os.remove(output_file)
				break
				
	# produce the initialization end message contains the configuration
	with topic_to_client.get_sync_producer() as producer:
		print 'Sending:'
		end_msg = auxiliary.generate_message(appliance_id=appliance_id, transaction_id=transaction_id, msg_type='initialization_end', value=output)
		producer.produce(end_msg)

	print 'Appliance initialization successful.'
