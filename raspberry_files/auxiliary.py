import time
import random
import json

def generate_message(appliance_id=None, msg_type=None, transaction_id=None, value=None):
	message = {}
	
	if msg_type:
		message['type'] = msg_type
		
	if value is not None:
		message['value'] = value
	
	if transaction_id:
		message['transaction_id'] = transaction_id
	
	if appliance_id:
		message['appliance_id'] = appliance_id
	
	message['timestamp'] = get_time()
	return json.dumps(message)


def read_config(config_path):
	with open(config_path) as config_file:    
		config = json.load(config_file)
		return config
		
def get_time():
	return str(int(time.time()))
