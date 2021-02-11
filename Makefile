REG := 100225593120.dkr.ecr.us-east-1.amazonaws.com
DOCKER_PULL_TAG := latest
DOCKER_BUILD_TAG := latest
ALLIANCE_RELEASE := 3.1.0

registry-docker-login:
ifneq ($(shell echo ${REG} | egrep "ecr\..+\.amazonaws\.com"),)
	@$(eval DOCKER_LOGIN_CMD=aws)
ifneq (${AWS_PROFILE},)
	@$(eval DOCKER_LOGIN_CMD=${DOCKER_LOGIN_CMD} --profile ${AWS_PROFILE})
endif
	@$(eval DOCKER_LOGIN_CMD=${DOCKER_LOGIN_CMD} ecr get-login-password | docker login -u AWS --password-stdin https://${REG})
	${DOCKER_LOGIN_CMD}
endif

build: pull
	docker build -t agrlocal/agr_file_generator:${DOCKER_BUILD_TAG} --build-arg REG=${REG} --build-arg DOCKER_PULL_TAG=${DOCKER_PULL_TAG} .

pull: registry-docker-login
	docker pull ${REG}/agr_base_linux_env:${DOCKER_PULL_TAG}

run: build
	DOCKER_BUILD_TAG=${DOCKER_BUILD_TAG} ALLIANCE_RELEASE=${ALLIANCE_RELEASE} docker-compose up agr_file_generator

test: run
	DOCKER_BUILD_TAG=${DOCKER_BUILD_TAG} ALLIANCE_RELEASE=${ALLIANCE_RELEASE} docker-compose up agr_file_generator_test
