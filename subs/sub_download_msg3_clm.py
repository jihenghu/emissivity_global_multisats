
import eumdac     ## pip3 istall eumdac
import datetime
import shutil
import requests
import os,sys

if len(sys.argv) < 3:
    print("ERROR! Non-adequate arguments for MSG Cloud Download Script !")
    exit()

yyyymmdd= sys.argv[1]
HHMM= sys.argv[2]
dirc= sys.argv[3]


# yyyymmdd="20221104"
# HHMM="1445"
# dirc='/data/jihenghu/data/MSG_CLM_IODC/'
## Download a single product

year=int(yyyymmdd[0:4])
month=int(yyyymmdd[4:6])
day=int(yyyymmdd[6:8])
hh=int(HHMM[0:2])
mm=int(HHMM[2:4])

dir_msg=dirc+'/'+yyyymmdd+'/'
if not os.path.exists(dirc):
    os.mkdir(dirc)
if not os.path.exists(dir_msg):
    os.mkdir(dir_msg)


# Insert your personal key and secret into the single quotes
consumer_key = 'LWSvREhJ9Gf9dfTURaFtP6aIV4Qa'
consumer_secret = '1fM914MR_biGBMmYtbfLPKihY7oa'

credentials = (consumer_key, consumer_secret)

token = eumdac.AccessToken(credentials)

try:
    pass# print(f"This token '{token}' expires {token.expiration}")
except requests.exceptions.HTTPError as error:
    print(f"Unexpected error: {error}")

datastore = eumdac.DataStore(token)

# get collection
try:    
    # https://data.eumetsat.int/product/EO:EUM:DAT:MSG:CLM-IODC 
    selected_collection = datastore.get_collection('EO:EUM:DAT:MSG:CLM') 
    print(f"{selected_collection} - {selected_collection.title}")
except eumdac.datastore.DataStoreError as error:
    print(f"Error related to the data store: '{error.msg}'")
except eumdac.collection.CollectionError as error:
    print(f"Error related to the collection: '{error.msg}'")
except requests.exceptions.ConnectionError as error:
    print(f"Error related to the connection: '{error.msg}'")
except requests.exceptions.RequestException as error:
    print(f"Unexpected error: {error}")

# Set sensing start and end time
start = datetime.datetime(year, month, day, hh, mm)
end = datetime.datetime(year, month, day, hh, mm)

# Retrieve datasets that match our filter
products = selected_collection.search(dtstart=start, dtend=end)

for product in products:
    try:
        with product.open(entry=str(product)+'.grb') as fsrc, \
                open(dir_msg+fsrc.name, mode='wb') as fdst:
            shutil.copyfileobj(fsrc, fdst)
            print(f'Download of file {fsrc.name} finished.')
    except eumdac.product.ProductError as error:
        print(f"Error related to the product '{product}' while trying to download it: '{error.msg}'")
    except requests.exceptions.ConnectionError as error:
        print(f"Error related to the connection: '{error.msg}'")
    except requests.exceptions.RequestException as error:
        print(f"Unexpected error: {error}")
