function [delta] = Flyby(velocity_in, h_p, time, bodies)

    GM = bodies.gm;
    RP = bodies.radius;
    state_vec_p = KM_solver(time, bodies);
    
    v_in = velocity_in - state_vec_p(4:6);
    delta = 2*asin((GM/(RP+h_p*RP)/(norm(v_in)^2+ GM/(RP+h_p*RP))));
    
end