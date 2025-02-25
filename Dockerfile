FROM mongo

WORKDIR /mnt/d/Rocksdb_LI/rocksdb
COPY ./ /src/rocksdb_li/

# EXPOSE 27017
# CMD ["mongod --storageEngine rocksdb"]
