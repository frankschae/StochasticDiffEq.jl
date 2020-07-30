@muladd function perform_step!(integrator,cache::EXEMConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,W,p = integrator

  K = uprev .+ dt .* integrator.f(uprev,p,t)

  if !is_diagonal_noise(integrator.sol.prob) || typeof(W.dW) <: Number
    noise = integrator.g(uprev,p,t)*W.dW
  else
    noise = integrator.g(uprev,p,t).*W.dW
  end
  u = K + noise

  integrator.u = u
end

@muladd function perform_step!(integrator,cache::EXEMCache,f=integrator.f)
  @unpack tmp,rtmp1,rtmp2 = cache
  @unpack t,dt,uprev,u,W,p = integrator
  integrator.f(rtmp1,uprev,p,t)

  @.. u = uprev + dt * rtmp1

  integrator.g(rtmp2,u_choice,p,t)

  if is_diagonal_noise(integrator.sol.prob)
    @.. rtmp2 *= W.dW
    @.. u = u + rtmp2
  else
    mul!(rtmp1,rtmp2,W.dW)
    @.. u = u + rtmp1
  end
end
