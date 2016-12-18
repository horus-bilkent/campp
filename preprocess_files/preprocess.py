#!/usr/bin/env python
import time
import json
import random
import sys
import jsonpickle
import os
import subprocess

from pykafka import KafkaClient, SslConfig

CONFIG_FILE='./preprocess_config.json'

def generate_initialization_message(appliance_id, value):
	message = {}
	message['appliance_id'] = appliance_id
	message['value'] = value
	message['timestamp'] = str(int(time.time()))
	return json.dumps(message)

def write_jpg(jpg_bytes, to_file):
	with open(to_file, 'w') as outfile:
		outfile.write(jsonpickle.decode(jpg_bytes))
	

if __name__ == "__main__":	
	with open(CONFIG_FILE) as config_file:    
		config = json.load(config_file)
		
	try:
		
		broker_connect = config['broker_connect']
		initialization_topic = config['initialization_topic']
		preprocess = config['preprocess']
	except KeyError, e:
		print 'Error in config file: ' + str(e)
		exit(1)
		
	# config = SslConfig(cafile=ssl_ca_cert, certfile=ssl_client_cert, keyfile=ssl_client_private_key_file, password=ssl_key_password)
	# client = KafkaClient(hosts=broker_connect, ssl_config=config)
	client = KafkaClient(hosts=broker_connect)
	
	topic_initialization = client.topics[bytes(initialization_topic)]
	
	consumer = topic_initialization.get_simple_consumer(fetch_message_max_bytes=50000000)
	for message_recv in consumer:
		if message_recv is not None:
			message_recv_json = json.loads(message_recv.value)
			if 'timestamp' in message_recv_json and int(time.time()) <= int(message_recv_json['timestamp']) + 15: 
				print 'Received.'
				image_1_bytes = message_recv_json['image_1']
				image_2_bytes = message_recv_json['image_2']
				image_1_name = './image_1_' + str(int(time.time())) + '.jpg'
				image_2_name = './image_2_' + str(int(time.time())) + '.jpg'
				output_file = './preprocess_output'
				
				write_jpg(image_1_bytes, image_1_name)
				write_jpg(image_2_bytes, image_2_name)
				ret_val = 0

				# ret_val = subprocess.call([preprocess, image_1_name, image_2_name, output_file])
				
				if ret_val != 0:
					print 'Initialization failed, cannot execute preprocessing'
					exit(1)
				
				reply_to = message_recv_json['reply_to']
				topic_to_client = client.topics[bytes(reply_to)]
				
				with open(output_file, 'r') as infile:
					output = infile.read()
					
					
				# os.remove(image_1_name)
				# os.remove(image_2_name)
				# os.remove(output_file)
				
				
				with topic_to_client.get_sync_producer() as producer:
					print 'Sending:'
					send_message = generate_initialization_message(message_recv_json['appliance_id'], output)
					producer.produce(send_message)
					break
			else:
				print 'Skipping old message'
					
