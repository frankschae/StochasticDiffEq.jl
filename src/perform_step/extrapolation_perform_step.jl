@muladd function perform_step!(integrator,cache::EXEMConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,W,p = integrator
  @unpack expcoef = cache

  order = integrator.alg.order
  expcoeforder = @view expcoef[:,order]
  u = zero(integrator.u)

  for i in 1:order
    dt_int = dt * i
    sqdt_int = sqrt(dt_int)

    _dW = map(x -> calc_twopoint_random(sqdt_int, x),  W.dW)

    K = uprev .+ dt_int .* integrator.f(uprev,p,t)

    if !is_diagonal_noise(integrator.sol.prob) || typeof(W.dW) <: Number
      noise = integrator.g(uprev,p,t)*_dW
    else
      noise = integrator.g(uprev,p,t).*_dW
    end
    #u += @.. expcoeforder[i]*integrator.alg.func(K + noise)
    u += @.. expcoeforder[i]*(K + noise)
  end

  integrator.u = u
end

@muladd function perform_step!(integrator,cache::EXEMCache,f=integrator.f)
  @unpack tmp,rtmp1,rtmp2 = cache
  @unpack t,dt,uprev,u,W,p = integrator
  integrator.f(rtmp1,uprev,p,t)

  @.. u = uprev + dt * rtmp1

  integrator.g(rtmp2,u,p,t)

  if is_diagonal_noise(integrator.sol.prob)
    @.. rtmp2 *= W.dW
    @.. u = u + rtmp2
  else
    mul!(rtmp1,rtmp2,W.dW)
    @.. u = u + rtmp1
  end
end
