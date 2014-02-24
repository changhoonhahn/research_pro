function get_qf, Mstar, z
    qf_z0 = -6.04 + 0.63*Mstar
    case 1 of
    (Mstar GE 9.5) AND (Mstar LT 10.0): alpha = -2.1
    (Mstar GE 10.0) AND (Mstar LT 10.5): alpha = -2.2
    (Mstar GE 10.5) AND (Mstar LT 11.0): alpha = -2.0
    (Mstar GE 11.0) AND (Mstar LT 11.5): alpha = -1.3
    else: print, "Mstar is out of range"
    endcase
    
    output = qf_z0*(1.0+z)^alpha 
    if (output LT 0.0) then output = 0.0
    return, output
end
