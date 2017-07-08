//
// Created by Vasiliy Ershov on 27/03/2017.
//

#include "normal_quality_model.hpp"

using namespace n_normal_model;

std::vector<double> NormalClusterModel::left_likelihoods_ = {
    -9.98, -4.95, -3.95, -3.5, -3, -2.5, -2.2, -2};
std::vector<double> NormalClusterModel::equal_likelihoods_ = {
    -0.001, -0.001, -0.019, -0.05, -0.07, -0.15, -0.2, -0.25};
std::vector<double> NormalClusterModel::right_likelihoods_ = {
    -5.99, -5.95, -5, -4.35, -3.8, -3, -2.8, -2.5};