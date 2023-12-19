FROM node:16-alpine

ARG HTTP_PORT=80
ARG HTTPS_PORT=443
ENV HTTP_PORT=$HTTP_PORT
ENV HTTPS_PORT=$HTTPS_PORTENV
ENV ARGS=''

WORKDIR /app


# Fixes permission issues when deployed to Openshift
# Similar fix added to Vault Dockerfile
ENV NPM_FOLDER=/.npm
RUN mkdir -p $NPM_FOLDER && \
    chgrp -R 0 $NPM_FOLDER && \
    chmod -R g=u $NPM_FOLDER

COPY server/package-lock.json server/package.json ./
RUN npm install --production
COPY server .
RUN cd ui && npm install && npm run build

EXPOSE ${HTTP_PORT} ${HTTPS_PORT}

CMD npm run start ${ARGS}
