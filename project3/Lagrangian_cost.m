function J=Lagrangian_cost(D,R,Q)

lamada=0.2*Q*Q;
J=D+lamada*R;
end