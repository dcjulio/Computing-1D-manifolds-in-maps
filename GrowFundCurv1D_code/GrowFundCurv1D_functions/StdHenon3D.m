classdef StdHenon3D
    methods     ( Static = true )

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%% MAPPING %%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function outpoints = mapping(inpoints,stab,opts)
            %manif.points: coordinates x,y,z
            %manif.stab: 'Umanifold' or 'Smanifold'
            points=inpoints;
            for i=1:opts.mapiter % time the map is applied

                %decompactify
                decomp_points=StdHenon3D.decompactify(points);

                % using map or inverse map depending on the stability of the manifold
                if isequal('Umanifold',stab)
                    map_points=StdHenon3D.ff(decomp_points,opts);

                elseif isequal('Smanifold',stab)
                    map_points=StdHenon3D.ff_inv(decomp_points,opts);

                else
                    fprintf('\nError!: Stability is not defined correctly. \nUse Smanifold for stable or Umanifold for unstable\n')
                    break; 
                end

                %compactify
                points=StdHenon3D.compactify(map_points); 

            end
            outpoints=points;
      
        end
   


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% DEFINITION OF THE MAP %%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        % > -------- map
        function map_points=ff(points,opts)

            a=opts.par.a;
            b=opts.par.b;
            xi=opts.par.xi;
            
            map_points=struct();

            % define the map
            map_points.x = points.y;
            map_points.y = a - points.y.^2 - b.*points.x;
            map_points.z = xi*points.z + points.y;

        end
        
        %----------------------------------------------
        %----------------------------------------------
        %----------------------------------------------

        % > -------- inverse map
        function map_points=ff_inv(points,opts)

            a=opts.par.a;
            b=opts.par.b;
            xi=opts.par.xi;
            
            map_points=struct();

            % define the inverse map
            map_points.x = (a - points.x.^2 - points.y)./b;
            map_points.y = points.x;
            map_points.z = (-points.x + points.z)/xi;

        end

        %----------------------------------------------
        %----------------------------------------------
        %----------------------------------------------




        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%% FIXED POINTS %%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % > -------- fixed points 
        function fixpinfo=fixpoints(opts)

            a=opts.par.a;
            b=opts.par.b;
            xi=opts.par.xi;

            fixpinfo=struct();

            % p minus fixed point
            fixpinfo.pmin.x=1/2*(-1-b-sqrt(4*a+(1+b)^2));
            fixpinfo.pmin.y=1/2*(-1-b-sqrt(4*a+(1+b)^2));
            fixpinfo.pmin.z=(1+b+sqrt(4*a+(1+b)^2))/(2*(-1+xi));

            % p plus fixed point
            fixpinfo.pplu.x=1/2*(-1-b+sqrt(4*a+(1+b)^2));
            fixpinfo.pplu.y=1/2*(-1-b+sqrt(4*a+(1+b)^2));
            fixpinfo.pplu.z=(1+b-sqrt(4*a+(1+b)^2))/(2*(-1+xi));

            

             % computing them in compactified coordinates
             fixp_names=fieldnames(fixpinfo);
             for k=1:numel(fixp_names) 
                 fixpinfo.(fixp_names{k})=StdHenon3D.compactify(fixpinfo.(fixp_names{k}));
             end
            

        end

         
        %------------ Compactification ---------------
        %----------------------------------------------
        %----------------------------------------------
        %----------------------------------------------
        
        
        % > -------- decompactify
        function decomp_points=decompactify(points)
            decomp_points=struct();

            r = (1 - points.x.^2 - points.y.^2);
            decomp_points.x = 2*points.x./r;
            decomp_points.y = 2*points.y./r;
            decomp_points.z = 2*points.z./(1 - points.z.^2);

        end
        
        %----------------------------------------------
        %----------------------------------------------
        %----------------------------------------------
    
        % > -------- compactify
        function comp_points=compactify(points)
            comp_points=struct();

            r = 1 + sqrt(1 + points.x.^2 + points.y.^2);
            comp_points.x = points.x./r;
            comp_points.y = points.y./r;
            comp_points.z = points.z./(1 + sqrt(1 + points.z.^2));

        end
        

        %----------------------------------------------
        %----------------------------------------------
        %----------------------------------------------
        % > -------- map and inverse map from compactification to compactification
        
        function inv_points=invmap_comp(points,opts)

            decomp=StdHenon3D.decompactify(points);
            inv=StdHenon3D.ff_inv(decomp,opts);
            inv_points=StdHenon3D.compactify(inv);

        end
        
        function inv_points=map_comp(points,opts)

            decomp=StdHenon3D.decompactify(points);
            inv=StdHenon3D.ff(decomp,opts);
            inv_points=StdHenon3D.compactify(inv);

        end

    end
end

