struct EXEMConstantCache <: StochasticDiffEqConstantCache end
@cache struct EXEMCache{uType,rateType,rateNoiseType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  rtmp1::rateType
  rtmp2::rateNoiseType
end

alg_cache(alg::EXEM,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}}) = EXEMConstantCache()

function alg_cache(alg::EXEM,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  tmp = zero(u); rtmp1 = zero(rate_prototype);
  if noise_rate_prototype !== nothing
    rtmp2 = zero(noise_rate_prototype)
  else
    rtmp2 = nothing
  end
  EXEMCache(u,uprev,tmp,rtmp1,rtmp2)
end
