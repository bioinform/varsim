#build image and push
version=$(git describe | sed 's/^v//')
branch=$(git rev-parse --abbrev-ref HEAD)
docker build --no-cache --build-arg branch_of_interest=${branch} -f Dockerfile -t carlocreator/varsim:${version} .
docker build --no-cache --build-arg branch_of_interest=${branch} -f Dockerfile -t carlocreator/varsim:latest .
docker push carlocreator/varsim:${version}
docker push carlocreator/varsim:latest
