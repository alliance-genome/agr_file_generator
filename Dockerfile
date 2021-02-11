ARG REG=100225593120.dkr.ecr.us-east-1.amazonaws.com
ARG DOCKER_IMAGE_TAG=latest

FROM ${REG}/agr_base_linux_env:${DOCKER_IMAGE_TAG}

WORKDIR /usr/src/app

ADD requirements.txt .

RUN pip3 install -r requirements.txt

RUN mkdir tmp

ADD . .

CMD ["python3", "-u", "src/app.py", "--all-filetypes", "--upload"]
