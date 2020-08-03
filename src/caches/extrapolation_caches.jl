struct EXEMConstantCache{SeqType} <: StochasticDiffEqConstantCache
  # hard-coded Romberg sequence
  R1::SeqType
  R2::SeqType
  R3::SeqType
  R4::SeqType
  R5::SeqType
end

function EXEMConstantCache(order, T::Type)

  R1 = convert(T, 1//2)
  R2 = convert(T, -1)
  R3 = convert(T, 2)
  R4 = convert(T, 342//491)
  R5 = convert(T, 342//491)

  EXEMConstantCache()
end


@cache struct EXEMCache{uType,rateType,rateNoiseType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  rtmp1::rateType
  rtmp2::rateNoiseType
end

function alg_cache(alg::EXEM,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}})
    if alg.order == 1
      R1::T
      R2::T
      R3::T
      R4::T
      R5::T
    elseif alg.order == 2

    elseif alg.order == 3

    elseif alg.order == 4

    elseif alg.order == 5

    end
    EXEMConstantCache(alg.order, real(uBottomEltypeNoUnits))
end

function alg_cache(alg::EXEM,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  tmp = zero(u); rtmp1 = zero(rate_prototype);
  if noise_rate_prototype !== nothing
    rtmp2 = zero(noise_rate_prototype)
  else
    rtmp2 = nothing
  end
  EXEMCache(u,uprev,tmp,rtmp1,rtmp2)
end
