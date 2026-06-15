cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(chill-plus VERSION 2.0.0 LANGUAGES CXX)

find_package(coretoolkit REQUIRED)

volt_add_plugin(chill-plus
    MAIN_SOURCE src/plugin.cpp
)
