// (C)  Copyright Alexander (Lex) Kemper akemper@ncsu.edu 2022


#pragma once

#include<tuple>
#include<optional>
#include<functional>

#include "cartan/algebra_types.h"
#include "cartan/involutions.h"

typedef std::function<bool(pauli_int)> involutionFunc;

std::tuple<algebra, algebra> get_algebra_by_commuting(
        algebra* p_algebra_g, 
        algebra* p_algebra_h = nullptr,
        const std::optional<involutionFunc>& involution = std::nullopt
        );

