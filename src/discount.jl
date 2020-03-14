
# Effective Rate Factor
ERF(::ContinuousCompounding, r::Float64, t::YearFraction) = iszero(t) ? 1.0 : exp(r*value(t))
ERF(::SimpleCompounding, r::Float64, t::YearFraction) = iszero(t) ? 1.0 : 1.0 + r*value(t)
ERF(::ExponentialCompounding, r::Float64, t::YearFraction) = iszero(t) ? 1.0 : (1.0+r)^value(t)
ERF(ct::CompoundingType, dcc::DayCountConvention, r::Float64, date_start::Date, date_end::Date) = ERF(ct, r, yearfraction(dcc, date_start, date_end))
ERF(curve::AbstractIRCurve, maturity::Date) = ERF(curve_get_compounding(curve), curve_get_daycount(curve), zero_rate(curve, maturity), curve_get_date(curve), maturity)
ERF(curve::AbstractIRCurve, maturity::YearFraction) = ERF(curve_get_compounding(curve), zero_rate(curve, maturity), maturity)

function ERF(curve::AbstractIRCurve, forward_date::Date, maturity::Date)
    @assert forward_date <= maturity "maturity $maturity should follow the forward date $forward_date."
    ERF(curve, maturity) / ERF(curve, forward_date)
end

function ERF(curve::AbstractIRCurve, forward_date::YearFraction, maturity::YearFraction)
    @assert forward_date <= maturity "maturity $maturity should follow the forward date $forward_date."
    ERF(curve, maturity) / ERF(curve, forward_date)
end

ERF_to_rate(::ContinuousCompounding, ERF::Float64, t::YearFraction) = log(ERF) / value(t)
ERF_to_rate(::SimpleCompounding, ERF::Float64, t::YearFraction) = (ERF - 1.0) / value(t)
ERF_to_rate(::ExponentialCompounding, ERF::Float64, t::YearFraction) = ERF^(1.0/value(t)) - 1.0
ERF_to_rate(curve::AbstractIRCurve, ERF::Float64, t::YearFraction) = ERF_to_rate(curve_get_compounding(curve), ERF, t)

discountfactor_to_rate(c::CompoundingType, _discountfactor_::Float64, t::YearFraction) = ERF_to_rate(c, 1.0 / _discountfactor_, t)

function discountfactor_to_rate(c::CompoundingType, _discountfactor_vec_::Vector{Float64}, t_vec::Vector{YearFraction})
    l = length(_discountfactor_vec_)
    @assert l == length(t_vec) "_discountfactor_vec_ and t_vec must have the same length. ($l != $(length(t_vec)))"

    result = Vector{Float64}(undef, l)
    for i in 1:l
        @inbounds result[i] = discountfactor_to_rate(c, _discountfactor_vec_[i], t_vec[i])
    end

    return result
end

# Effective Rate = [Effective Rate Factor] - 1
ER(c::CompoundingType, r::Float64, t::YearFraction) = ERF(c, r, t) - 1.0
ER(ct::CompoundingType, dcc::DayCountConvention, r::Float64, date_start::Date, date_end::Date) = ER(ct, r, yearfraction(dcc, date_start, date_end))
ER(curve::AbstractIRCurve, maturity::Date) = ERF(curve, maturity) - 1.0
ER(curve::AbstractIRCurve, forward_date::Date, maturity::Date) =  ERF(curve, forward_date, maturity) - 1.0

# [Discount Factor] = 1 / [Effective Rate Factor]
discountfactor(ct::CompoundingType, r::Float64, t::YearFraction) = 1.0 / ERF(ct, r, t)
discountfactor(ct::CompoundingType, dcc::DayCountConvention, r::Float64, date_start::Date, date_end::Date) = discountfactor(ct, r, yearfraction(dcc, date_start, date_end))

# maturity can be either a date or an year fraction
function discountfactor(curve::AbstractIRCurve, maturity::T) where {T<:Union{Date, YearFraction}}
    1.0 / ERF(curve, maturity)
end

function discountfactor(curve::AbstractIRCurve, forward_date::T1, maturity::T2) where {T1<:Union{Date, YearFraction}, T2<:Union{Date, YearFraction}}
    1.0 / ERF(curve, forward_date, maturity)
end

# Optimized vector functions for `ERF` and `discountfactor` functions
for fun in (:ERF, :discountfactor)
    @eval begin
        function ($fun)(curve::AbstractIRCurve, maturity_vec::Vector{Date})
            len = length(maturity_vec)
            _zero_rate_vec_ = zero_rate(curve, maturity_vec)
            result = Vector{Float64}(undef, len)
            for i = 1:len
                @inbounds result[i] = $(fun)(curve_get_compounding(curve), curve_get_daycount(curve), _zero_rate_vec_[i], curve_get_date(curve), maturity_vec[i])
            end
            return result
        end
    end
end
