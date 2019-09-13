#!/bin/bash
# docker pull quay.io/pypa/manylinux1_x86_64
# # or docker pull quay.io/pypa/manylinux1_i686
# docker run -it -v (pwd):/io quay.io/pypa/manylinux1_x86_64

for PYBIN in /opt/python/*3[5-7]*/bin; do
    "${PYBIN}/pip" wheel /io/ -w wheelhouse/
done

for whl in wheelhouse/*.whl; do
    auditwheel repair "$whl" -w /io/wheelhouse/
done
