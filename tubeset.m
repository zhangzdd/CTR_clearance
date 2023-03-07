function tset = tubeset(tubes,ctrl)
    tset.ctrl = ctrl;% Control variables {[base_totation_angle,base_translation_length],...}
    tset.tube = tubes;
    tset.cnt = size(tubes,2);
end