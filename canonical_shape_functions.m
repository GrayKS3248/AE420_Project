    % Volume calculation
    volume = det([1 P1'; 1 P2'; 1 P3'; 1 P4']) / 6.0;
    
    % Shape function 1 constants calculation
    a1 = det([P2'; P3'; P4']);
    b1 = -det([1 P2(2:3)'; 1 P3(2:3)'; 1 P4(2:3)']);
    c1 = det([1 P2(1) P2(3); 1 P3(1) P3(3); 1 P4(1) P4(3)]);
    d1 = -det([1 P2(1:2)'; 1 P3(1:2)'; 1 P4(1:2)']);
    
    % Shape function 2 constants calculation
    a2 = -det([P1'; P3'; P4']);
    b2 = det([1 P1(2:3)'; 1 P3(2:3)'; 1 P4(2:3)']);
    c2 = -det([1 P1(1) P1(3); 1 P3(1) P3(3); 1 P4(1) P4(3)]);
    d2 = det([1 P1(1:2)'; 1 P3(1:2)'; 1 P4(1:2)']);
    
    % Shape function 3 constants calculation
    a3 = det([P1'; P2'; P4']);
    b3 = -det([1 P1(2:3)'; 1 P2(2:3)'; 1 P4(2:3)']);
    c3 = det([1 P1(1) P1(3); 1 P2(1) P2(3); 1 P4(1) P4(3)]);
    d3 = -det([1 P1(1:2)'; 1 P2(1:2)'; 1 P4(1:2)']);
    
    % Shape function 4 constants calculation
    a4 = -det([P1'; P2'; P3']);
    b4 = det([1 P1(2:3)'; 1 P2(2:3)'; 1 P3(2:3)']);
    c4 = -det([1 P1(1) P1(3); 1 P2(1) P2(3); 1 P3(1) P3(3)]);
    d4 = det([1 P1(1:2)'; 1 P2(1:2)'; 1 P3(1:2)']);

    cdN = (1.0/(6.0*volume)) * [b1  0   0   b2  0   0   b3  0   0   b4  0   0;
                                0   c1  0   0   c2  0   0   c3  0   0   c4  0; 
                                0   0   d1  0   0   d2  0   0   d3  0   0   d4;
                                c1  b1  0   c2  b2  0   c3  b3  0   c4  b4  0;
                                0   d1  c1  0   d2  c2  0   d3  c3  0   d4  c4;
                                d1  0   b1  d2  0   b2  d3  0   b3  d4  0   b4];