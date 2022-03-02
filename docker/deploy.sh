#!/bin/bash

echo "./deploy.sh $*" > redeploy.sh
chmod +x redeploy.sh

if [ -z "$(docker network ls -qf name=data-table-network)" ]; then
    echo "creating network"
    docker network create --driver bridge data-table-network
fi

existing=$(docker ps -aqf name='data-table(-server)?$')
if [ -n "$existing" ]; then
    echo "removing existing container(s)"
    docker rm -f $existing
fi

# NGINX and SERVER_ARGS used for jwilder/nginx-proxy config
NGINX=0
PORT_SET=0
SERVER_ARGS=()
DEFAULT_PORT=80
SERVER_PORT=80
PORT=
ARGS=$*

while [ -n "$1" ]; do
    case $1 in
        --https )
            if [ -z "$PORT" ]; then
                PORT=443
            fi
            SERVER_PORT=443
            SERVER_ARGS+=(
                --env CERT_NAME=default
                --env VIRTUAL_PORT=443
                --env VIRTUAL_PROTO=https
            )
            ;;
        --nginx )
            NGINX=1
            ;;
        -u | --url )
            SERVER_ARGS+=(--env VIRTUAL_HOST=$2)
            shift
            ;;
        -w | --web-port )
            PORT=$2
            PORT_SET=1
            shift
            ;;
    esac
    shift
done

if [ -z "$PORT" ]; then
    PORT=$DEFAULT_PORT
elif [ $NGINX -eq 1 ] && [ $PORT_SET -eq 1 ]; then
    echo "Error: --nginx and -w/--web-port are mutually exclusive"
    exit 1
fi

# bind port if not using nginx proxy
if [ $NGINX -eq 0 ]; then
    SERVER_ARGS+=(-p $PORT:$SERVER_PORT)
fi

docker run -d \
--name data-table \
--restart unless-stopped \
--network data-table-network \
--env no_proxy=data-table-server \
--env NO_PROXY=data-table-server \
-e ARGS="$ARGS" \
data-table

docker run -d \
--name data-table-server \
--restart unless-stopped \
--network data-table-network \
${SERVER_ARGS[@]} \
data-table-server
