#!/bin/bash
set -euxo pipefail

# Please execute this script from within the root apex directory.
if [[ ! -f ".dockerignore" ]]; then
  echo "Must execute from apex root directory."
  exit 1
fi

APEX_VERSION=`git describe --tags --abbrev=11 | sed 's/^v//' | sed 's/-g/-/'`

# Build docker image (which builds a linux static executable)
bin/docker_build_image.sh

# Copy executable out of docker container
docker run --rm apex:latest cat /home/apex/build/apex > _apex

# Create directory for distribution
RELEASE_DIR="apex-${APEX_VERSION}-x86_64-unknown-linux-gnu"
mkdir -p "dist/${RELEASE_DIR}"

# Move executable
mv _apex "dist/${RELEASE_DIR}/apex"

# Copy other necessary files
# TODO: Corbin should designate a license; the LICENSE file should be included in the tarball here
cp README.md "dist/${RELEASE_DIR}"

# Create tarball
cd dist
tar zcf ${RELEASE_DIR}.tar.gz ${RELEASE_DIR}
