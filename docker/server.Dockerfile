FROM node:16-alpine

ARG HTTP_PORT=80
ARG HTTPS_PORT=443
ENV HTTP_PORT=$HTTP_PORT
ENV HTTPS_PORT=$HTTPS_PORTENV
ENV ARGS=''

WORKDIR /app

COPY server/package-lock.json server/package.json ./
RUN npm install --production
COPY server .
RUN cd ui && npm install && npm run build

EXPOSE ${HTTP_PORT} ${HTTPS_PORT}

CMD npm run start ${ARGS}
