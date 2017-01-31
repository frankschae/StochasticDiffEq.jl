# Use Recursion to find the first callback for type-stability

# Base Case: Only one callback
function find_first_continuous_callback(integrator,callback::ContinuousCallback)
  (find_callback_time(integrator,callback)...,1,1)
end

# Starting Case: Compute on the first callback
function find_first_continuous_callback(integrator,callback::ContinuousCallback,args...)
  find_first_continuous_callback(integrator,find_callback_time(integrator,callback)...,1,1,args...)
end

function find_first_continuous_callback(integrator,tmin::Number,upcrossing::Float64,idx::Int,counter::Int,callback2)
  counter += 1 # counter is idx for callback2.
  tmin2,upcrossing2 = find_callback_time(integrator,callback2)
  if tmin < tmin2
    return tmin,upcrossing,idx,counter
  else
    return tmin2,upcrossing,counter,counter
  end
end

function find_first_continuous_callback(integrator,tmin::Number,upcrossing::Float64,idx::Int,counter::Int,callback2,args...)
  find_first_continuous_callback(integrator,find_first_continuous_callback(integrator,tmin,upcrossing,idx,counter,callback2)...,args...)
end

@inline function determine_event_occurance(integrator,callback)
  event_occurred = false
  Θs = linspace(typeof(integrator.t)(0),typeof(integrator.t)(1),callback.interp_points)
  interp_index = 0
  # Check if the event occured
  previous_condition = callback.condition(integrator.tprev,integrator.uprev,integrator)
  if isapprox(previous_condition,0,rtol=callback.reltol,atol=callback.abstol)
    prev_sign = 0.0
  else
    prev_sign = sign(previous_condition)
  end
  prev_sign_index = 1
  if ((prev_sign<0 && !(typeof(callback.affect!)<:Void)) || (prev_sign>0 && !(typeof(callback.affect_neg!)<:Void))) && prev_sign*sign(callback.condition(integrator.tprev+integrator.dt,integrator.u,integrator))<0
    event_occurred = true
    interp_index = callback.interp_points
  elseif callback.interp_points!=0 # Use the interpolants for safety checking
    for i in 2:length(Θs)-1
      new_sign = callback.condition(integrator.tprev+integrator.dt*Θs[i],sde_interpolant(Θs[i],integrator),integrator)
      if prev_sign == 0
        prev_sign = new_sign
        prev_sign_index = i
      end
      if ((prev_sign<0 && !(typeof(callback.affect!)<:Void)) || (prev_sign>0 && !(typeof(callback.affect_neg!)<:Void))) && prev_sign*new_sign<0
        event_occurred = true
        interp_index = i
        break
      end
    end
  end
  event_occurred,interp_index,Θs,prev_sign,prev_sign_index
end

function find_callback_time(integrator,callback)
  event_occurred,interp_index,Θs,prev_sign,prev_sign_index = determine_event_occurance(integrator,callback)
  if event_occurred
    if typeof(callback.condition) <: Void
      new_t = zero(typeof(integrator.t))
    else
      if callback.interp_points!=0
        top_Θ = Θs[interp_index] # Top at the smallest
      else
        top_Θ = typeof(integrator.t)(1)
      end
      if callback.rootfind
        find_zero = (Θ) -> begin
          callback.condition(integrator.tprev+Θ*integrator.dt,sde_interpolant(Θ,integrator),integrator)
        end
        Θ = prevfloat(prevfloat(fzero(find_zero,typeof(integrator.t)(0),top_Θ)))
        # 2 prevfloat guerentees that the new time is either 1 or 2 floating point
        # numbers just before the event, but not after. If there's a barrier
        # which is never supposed to be crossed, then this will ensure that
        # The item never leaves the domain. Otherwise Roots.jl can return
        # a float which is slightly after, making it out of the domain, causing
        # havoc.
        new_t = integrator.dt*Θ
      else
        # If no solve and no interpolants, just use endpoint
        new_t = integrator.dt
      end
    end
  else
    new_t = zero(typeof(integrator.t))
  end
  new_t,prev_sign
end

function apply_callback!(integrator,callback::ContinuousCallback,cb_time,prev_sign)
  if cb_time != zero(typeof(integrator.t))
    change_t_via_interpolation!(integrator,integrator.tprev+cb_time)
  end

  if callback.save_positions[1]
    update_running_noise!(integrator)
    savevalues!(integrator)
  end

  integrator.u_modified = true

  if prev_sign < 0
    if typeof(callback.affect!) <: Void
      integrator.u_modified = false
    else
      callback.affect!(integrator)
    end
  elseif prev_sign > 0
    if typeof(callback.affect_neg!) <: Void
      integrator.u_modified = false
    else
      callback.affect_neg!(integrator)
    end
  end

  if integrator.u_modified
    #reeval_internals_due_to_modification!(integrator)
    if callback.save_positions[2]
      if !callback.save_positions[1]
        update_running_noise!(integrator)
      end
      savevalues!(integrator)
    end
    return true
  end
  false
end

#Base Case: Just one
function apply_discrete_callback!(integrator::SDEIntegrator,callback::DiscreteCallback)
  if callback.save_positions[1]
    savevalues!(integrator)
  end

  integrator.u_modified = true
  if callback.condition(integrator.tprev+integrator.dt,integrator.u,integrator)
    callback.affect!(integrator)
    if callback.save_positions[2]
      savevalues!(integrator)
    end
  end
  integrator.u_modified
end

#Starting: Get bool from first and do next
function apply_discrete_callback!(integrator::SDEIntegrator,callback::DiscreteCallback,args...)
  apply_discrete_callback!(integrator,apply_discrete_callback!(integrator,callback),args...)
end

function apply_discrete_callback!(integrator::SDEIntegrator,discrete_modified::Bool,callback::DiscreteCallback,args...)
  bool = apply_discrete_callback!(integrator,apply_discrete_callback!(integrator,callback),args...)
  discrete_modified || bool
end

function apply_discrete_callback!(integrator::SDEIntegrator,discrete_modified::Bool,callback::DiscreteCallback)
  bool = apply_discrete_callback!(integrator,callback)
  discrete_modified || bool
end

resize!(integrator::SDEIntegrator,i::Int) = resize!(integrator,integrator.cache,i)
function resize!(integrator::SDEIntegrator,cache,i)
  prev_len = length(integrator.u)
  for c in full_cache(integrator)
    resize!(c,i)
  end
  for c in integrator.S₁
    resize!(c[2],i)
    resize!(c[3],i)
    if i > prev_len # fill in rands
      resize_noise_caches!(integrator,c,c[1],prev_len:i)
    end
  end
  for c in integrator.S₂
    resize!(c[2],i)
    resize!(c[3],i)
    if i > prev_len # fill in rands
      resize_noise_caches!(integrator,c,c[1],prev_len:i)
    end
  end
  resize!(integrator.ΔW,i)
  resize!(integrator.ΔZ,i)
  resize!(integrator.ΔWtilde,i)
  resize!(integrator.ΔZtilde,i)
  resize!(integrator.ΔWtmp,i)
  resize!(integrator.ΔZtmp,i)
  resize!(integrator.W,i)
  resize!(integrator.Z,i)
  if i > prev_len # fill in rands
    fill!(@view(integrator.W[prev_len:i]),zero(eltype(integrator.u)))
    fill!(@view(integrator.Z[prev_len:i]),zero(eltype(integrator.u)))
  end
end

function deleteat!(integrator::SDEIntegrator,i::Int)
  for c in full_cache(integrator)
    deleteat!(c,i)
  end
end

function terminate!(integrator::SDEIntegrator)
  integrator.opts.tstops.valtree = typeof(integrator.opts.tstops.valtree)()
end