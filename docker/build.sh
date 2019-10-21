#build image and push
version=$(git describe | sed 's/^v//')
branch=$(git rev-parse --abbrev-ref HEAD)
docker build --no-cache --build-arg branch_of_interest=${branch} -f Dockerfile -t rssbred/varsim:${version} .
docker build --no-cache --build-arg branch_of_interest=${branch} -f Dockerfile -t rssbred/varsim:latest .
docker push rssbred/varsim:${version}
docker push rssbred/varsim:latest
