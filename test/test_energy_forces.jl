# Tests energy and forces


for p = 3:9
   h = 0.1^p
   hVh = zeros(hV)
   for n = 1:length(matR)
      matR[n] += h
      r = norm.(R)
      dVh = mat(evaluate_d(eam, r, R))[:]
      hVh[:, n] = (dVh - dV) / h
      matR[n] -= h
   end
   @printf("%1.1e | %4.2e \n", h, vecnorm(hVh - hV, Inf))
end
