#!/bin/bash

VELOCITY_SET="D2Q9"
SIM_ID="000"

SRC="Shan-Chen.cpp"
EXE_NAME="${VELOCITY_SET}_${SIM_ID}"
EXE_DIR="../bin/${VELOCITY_SET}"
EXE_PATH="${EXE_DIR}/${EXE_NAME}"

CXX=g++
CXXFLAGS="-O3 -std=c++11"

mkdir -p "$EXE_DIR"

echo "Compilando $SRC -> $EXE_PATH..."
$CXX $CXXFLAGS "$SRC" -o "$EXE_PATH"

if [ $? -eq 0 ]; then
    echo "Compilação concluída com sucesso: $EXE_PATH"
else
    echo "Erro na compilação."
    exit 1
fi
