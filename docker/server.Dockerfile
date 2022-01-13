FROM node:10-alpine

ENV ARGS=''
WORKDIR /app

COPY server/package-lock.json server/package.json ./
RUN npm install --production
COPY server .

EXPOSE 80

CMD npm run start ${ARGS}
