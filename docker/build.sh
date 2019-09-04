#build image and push
version=$(git describe | sed 's/^v//')
docker build --no-cache -f Dockerfile -t carlocreator/varsim:${version} .
docker build --no-cache -f Dockerfile -t carlocreator/varsim:latest .
docker push carlocreator/varsim:${version}
docker push carlocreator/varsim:latest
