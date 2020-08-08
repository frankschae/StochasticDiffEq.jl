struct EXEMConstantCache{MType} <: StochasticDiffEqConstantCache
  # hard-coded Romberg sequence as a static matrix
  expcoef::MType
end

function EXEMConstantCache(max_order,T::Type)
  expcoef = SMatrix{max_order,max_order,T,max_order*max_order}(1, 0, 0, 0, 0,
                              2, -1, 0, 0, 0,
                              8//3, -6//3, 1//3, 0,0,
                              64//21,-56//21, 14//21, -1//21, 0,
                              1024//315,-960//315, 280//315,-30//315, 1//315)
  EXEMConstantCache(expcoef)
end


@cache struct EXEMCache{uType,randType,MType,rateType,rateNoiseType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  _dW::randType
  expcoef::MType
  rtmp1::rateType
  rtmp2::rateNoiseType
end

function alg_cache(alg::EXEM,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}})
    alg.order > alg.max_order && error("The maximum extrapolation order is: ", alg.max_order)
    EXEMConstantCache(alg.max_order,real(uBottomEltypeNoUnits))
end

function alg_cache(alg::EXEM,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})

  alg.order > alg.max_order && error("The maximum extrapolation order is: ", alg.max_order)

  if typeof(ΔW) <: Union{SArray,Number}
    _dW = copy(ΔW)
  else
    _dW = zero(ΔW)
  end

  tmp = zero(u); rtmp1 = zero(rate_prototype);
  if noise_rate_prototype !== nothing
    rtmp2 = zero(noise_rate_prototype)
  else
    rtmp2 = nothing
  end

  expcoef = EXEMConstantCache(alg.max_order,real(uBottomEltypeNoUnits))

  EXEMCache(u,uprev,tmp,_dW,expcoef,rtmp1,rtmp2)
end
