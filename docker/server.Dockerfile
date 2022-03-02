FROM node:16-alpine

ENV ARGS=''
WORKDIR /app

COPY server/package-lock.json server/package.json ./
RUN npm install --production
COPY server .
RUN cd ui && npm install && npm run build

EXPOSE 80 443

CMD npm run start ${ARGS}
