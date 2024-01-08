library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

-- Uncomment the following library declaration if using
-- arithmetic functions with Signed or Unsigned values
--use IEEE.NUMERIC_STD.ALL;

-- Uncomment the following library declaration if instantiating
-- any Xilinx leaf cells in this code.
library UNISIM;
use UNISIM.VComponents.all;

entity optimised_accurate is
generic (N:integer; M : integer); --N:width, M: heigth
port(
    a : in std_logic_vector(N-1 downto 0); --width of the multiplier
    b : in std_logic_vector(M-1 downto 0); --heigth of the multiplier
    d : in std_logic_vector(N downto 0); --has always to be zero, only necessary since generic generate
    p : out std_logic_vector(M+N-1 downto 0)
);
end optimised_accurate;

architecture Behavioral of optimised_accurate is
component LUT6_2
generic (INIT : bit_vector(63 downto 0) := X"0000000000000000");
port(
    O5 : out std_logic;
    O6 : out std_logic;
    I0 : in std_logic;
    I1 : in std_logic;
    I2 : in std_logic;
    I3 : in std_logic;
    I4 : in std_logic;
    I5 : in std_logic
);
end component;

component CARRY4
port(
    CO : out std_logic_vector(3 downto 0); 
    O : out std_logic_vector(3 downto 0);
    CI : in std_logic;
    CYINIT : in std_logic;
    DI : in std_logic_vector(3 downto 0);
    S : in std_logic_vector(3 downto 0)
);
end component;

--constant N : integer := 16;

--lut_o(stage)(lut position in stage), output of all luts, summarized per stage
type tLut_o is array(M/2-1 downto 0) of std_logic_vector(N+8 downto 0); --N+4
signal lut_o : tLut_o;
--carry(stage)(bit position in stage), carry word per stage including the carry in at position carry(stage)(0)
type tCarry is array(M/2-1 downto 0) of std_logic_vector(N+8 downto 0); --N+4
signal carry : tCarry;
--ccout(stage)(bit position in stage), carry chain output word per stage
type tOut is array(M/2-1 downto 0) of std_logic_vector(N+8 downto 0); --N+4
signal ccout : tOut;
--din(stage)(lut input), signal that is propagated to LUT AND carry chain, therefore additional word for generic carry chain generation
type tDin is array(M/2-1 downto 0) of std_logic_vector(N+7 downto 0); --N+3
signal din : tDin;

begin
--generic LUT A generation, special cases:
--row zero has no LUT at position zero
--row zero b(m-1) input is static '0'
--every other row LUT at position zero has all '0' as input
--every row LUT at position one input for a(n-1) is static '0'
--every row LUT at position N input for a(n) and a(n-1) is static a(N-1) 
GEN_LUTA :
for stage in 0 to M/2-1 generate
    LUTA0 :
    if stage = 0 generate
    GEN_LUTA0 : 
    for nr in 0 to N+1 generate
        LUTA00 : 
        if nr = 0 generate
            LUTA00X : LUT6_2
            generic map(
                INIT => X"0000007000000088"
            )
            port map(
                I0 => a(0),
                I1 => b(2*stage),
                I2 => b(2*stage+1),
                I3 => '0',
                I4 => '0',
                I5 => '1',
                O5 => ccout(stage)(nr),
                O6 => carry(stage)(nr)
            );
            end generate LUTA00;          
        LUTA0Y :
        if nr > 0 and nr < N generate
            LUTA0YX : LUT6_2
            generic map(
                INIT => X"AAA9655AA5596AAA"
            )
            port map(
                I0 => '0', --d(nr),
                I1 => a(nr-1),
                I2 => '0',
                I3 => b(2*stage),
                I4 => b(2*stage+1),
                I5 => a(nr),
                O6 => lut_o(stage)(nr-1)
             );
        end generate LUTA0Y;
        LUTAN : 
        if nr = N generate
            LUTB0N1 : LUT6_2
            generic map(
                INIT => X"00000000A999666A"
            )
            port map(          
                I0 => '0', --d(nr),
                I1 => a(nr-1),
                I2 => '0',
                I3 => b(2*stage),
                I4 => b(2*stage+1),
                I5 => '0',
                O6 => lut_o(stage)(nr-1)
            );
        end generate LUTAN;       
        LUTB0N1 : 
        if nr = N+1 generate
            LUTB0N1 : LUT6_2
            generic map(
                INIT => X"0000000056669995"
            )
            port map(
                I0 => '1',
                I1 => a(N-1),
                I2 => '0',
                I3 => b(2*stage),
                I4 => b(2*stage+1),
                I5 => '0',
                O6 => lut_o(stage)(nr-1)
            );
        end generate LUTB0N1;       
    end generate GEN_LUTA0;
    end generate LUTA0;
    LUTAY :
    if stage > 0 and stage < M/2 generate
    GEN_LUTAY : 
    for nr in 0 to N+1 generate
        LUTAY0 :
        if nr = 0 generate
            LUTAY0X : LUT6_2
            generic map(
                INIT => X"0BBF0880A66AA66A"
            )
            port map(
                I0 => ccout(stage-1)(nr+2),
                I1 => a(0),
                I2 => b(2*stage-1),
                I3 => b(2*stage),
                I4 => b(2*stage+1),
                I5 => '1',
                O5 => ccout(stage)(nr),
                O6 => carry(stage)(nr)
            );
        end generate LUTAY0;
        LUTAY1 :
        if nr > 0 and nr < N generate
            LUTAYYX : LUT6_2
            generic map(
                INIT => X"AAA9655AA5596AAA"
            )
            port map(
                I0 => ccout(stage-1)(nr+2),
                I1 => a(nr-1),
                I2 => b(2*stage-1),
                I3 => b(2*stage),
                I4 => b(2*stage+1),
                I5 => a(nr),
                O6 => lut_o(stage)(nr-1)
            );
        end generate LUTAY1;              
         LUTAYN :
         if  nr = N generate
            LUTAYNX : LUT6_2
            generic map(
                INIT => X"5666999555555555"
            )
            port map(
                I0 => carry(stage-1)(nr+1),
                I1 => a(nr-1),
                I2 => b(2*stage-1),
                I3 => b(2*stage),
                I4 => b(2*stage+1),
                I5 => '1',
                O5 => lut_o(stage)(nr+1),
                O6 => lut_o(stage)(nr-1)
            );
        end generate LUTAYN;
        LUTBYN1 : 
        if nr = N+1 generate
            LUTBYN1X : LUT6_2
            generic map(
                INIT => X"0000000056669995"
            )
            port map(
                I0 => carry(stage-1)(nr),
                I1 => a(nr-2),
                I2 => b(2*stage-1),
                I3 => b(2*stage),
                I4 => b(2*stage+1),
                I5 => '0',
                O6 => lut_o(stage)(nr-1)
            );
        end generate LUTBYN1;                           
     end generate GEN_LUTAY;
     end generate LUTAY;
end generate GEN_LUTA;

GEN_CARRY :
for stage in 0 to M/2-1 generate
    GEN_CARRYY : 
    for nr in 0 to N/4+1 generate --N/4 generate
    CARRYYX : CARRY4
    port map(
        CO => carry(stage) (nr*4+4 downto nr*4+1),
        O => ccout(stage) (nr*4+4 downto nr*4+1),
        CI => carry(stage) (nr*4),
        CYINIT => '0',
        DI => din(stage) (nr*4+3 downto nr*4),
        S => lut_o(stage) (nr*4+3 downto nr*4)                
    );
    end generate GEN_CARRYY;    
end generate GEN_CARRY;    







------------------------------------------------------------------------------------------
--Input & output configuration for 4x4
------------------------------------------------------------------------------------------
--din(0)(N downto 0) <= ('1' & d(N downto 1));
--din(1)(N downto 0) <= (carry(0)(N+1) & lut_o(1)(N+1) & ccout(0)(N+1 downto 3));

--p <= (ccout(1)(N+1 downto 0) & ccout(0)(1 downto 0));
------------------------------------------------------------------------------------------

------------------------------------------------------------------------------------------
--Input & output configuration for 6x6
------------------------------------------------------------------------------------------
--din(0)(N downto 0) <= ('1' & ( N downto 1 => '0'));--.d(N downto 1));
--din(1)(N downto 0) <= (carry(0)(N+1) & lut_o(1)(N+1) & ccout(0)(N+1 downto 3));
--din(2)(N downto 0) <= (carry(1)(N+1) & lut_o(2)(N+1) & ccout(1)(N+1 downto 3));

--p <= (ccout(2)(N+1 downto 0) & ccout(1)(1 downto 0) & ccout(0)(1 downto 0));
------------------------------------------------------------------------------------------

------------------------------------------------------------------------------------------
--Input & output configuration for 8x8
------------------------------------------------------------------------------------------
din(0)(N downto 0) <= ('1' & ( N downto 1 => '0'));--.d(N downto 1));
din(1)(N downto 0) <= (carry(0)(N+1) & lut_o(1)(N+1) & ccout(0)(N+1 downto 3));
din(2)(N downto 0) <= (carry(1)(N+1) & lut_o(2)(N+1) & ccout(1)(N+1 downto 3));
din(3)(N downto 0) <= (carry(2)(N+1) & lut_o(3)(N+1) & ccout(2)(N+1 downto 3));

p <= (ccout(3)(N+1 downto 0) & ccout(2)(1 downto 0) & ccout(1)(1 downto 0) & ccout(0)(1 downto 0));
------------------------------------------------------------------------------------------

------------------------------------------------------------------------------------------
--Input & output configuration for 8x4
------------------------------------------------------------------------------------------
--din(0)(N downto 0) <= ('1' & ( N downto 1 => '0'));--.d(N downto 1));
--din(1)(N downto 0) <= (carry(0)(N+1) & lut_o(1)(N+1) & ccout(0)(N+1 downto 3));

--p <= (ccout(1)(N+1 downto 0) & ccout(0)(1 downto 0));
------------------------------------------------------------------------------------------

------------------------------------------------------------------------------------------
--Input & output configuration for 10x10
------------------------------------------------------------------------------------------
--din(0)(N downto 0) <= ('1' & ( N downto 1 => '0'));--.d(N downto 1));
--din(1)(N downto 0) <= (carry(0)(N+1) & lut_o(1)(N+1) & ccout(0)(N+1 downto 3));
--din(2)(N downto 0) <= (carry(1)(N+1) & lut_o(2)(N+1) & ccout(1)(N+1 downto 3));
--din(3)(N downto 0) <= (carry(2)(N+1) & lut_o(3)(N+1) & ccout(2)(N+1 downto 3));
--din(4)(N downto 0) <= (carry(3)(N+1) & lut_o(4)(N+1) & ccout(3)(N+1 downto 3));

--p <= (ccout(4)(N+1 downto 0) & ccout(3)(1 downto 0) & ccout(2)(1 downto 0) & ccout(1)(1 downto 0) & ccout(0)(1 downto 0));
------------------------------------------------------------------------------------------

------------------------------------------------------------------------------------------
--Input & output configuration for 16x16
------------------------------------------------------------------------------------------
--din(0)(N downto 0) <= ('1' & ( N downto 1 => '0'));--.d(N downto 1));
--din(1)(N downto 0) <= (carry(0)(N+1) & lut_o(1)(N+1) & ccout(0)(N+1 downto 3));
--din(2)(N downto 0) <= (carry(1)(N+1) & lut_o(2)(N+1) & ccout(1)(N+1 downto 3));
--din(3)(N downto 0) <= (carry(2)(N+1) & lut_o(3)(N+1) & ccout(2)(N+1 downto 3));
--din(4)(N downto 0) <= (carry(3)(N+1) & lut_o(4)(N+1) & ccout(3)(N+1 downto 3));
--din(5)(N downto 0) <= (carry(4)(N+1) & lut_o(5)(N+1) & ccout(4)(N+1 downto 3));
--din(6)(N downto 0) <= (carry(5)(N+1) & lut_o(6)(N+1) & ccout(5)(N+1 downto 3));
--din(7)(N downto 0) <= (carry(6)(N+1) & lut_o(7)(N+1) & ccout(6)(N+1 downto 3));

--p <= (ccout(7)(N+1 downto 0) & ccout(6)(1 downto 0) & ccout(5)(1 downto 0) & ccout(4)(1 downto 0) & 
--    ccout(3)(1 downto 0) & ccout(2)(1 downto 0) & ccout(1)(1 downto 0) & ccout(0)(1 downto 0));
------------------------------------------------------------------------------------------

------------------------------------------------------------------------------------------
--Input & output configuration for 32x32
------------------------------------------------------------------------------------------
--din(0)(N downto 0) <= ('1' & ( N downto 1 => '0'));--.d(N downto 1));
--din(1)(N downto 0) <= (carry(0)(N+1) & lut_o(1)(N+1) & ccout(0)(N+1 downto 3));
--din(2)(N downto 0) <= (carry(1)(N+1) & lut_o(2)(N+1) & ccout(1)(N+1 downto 3));
--din(3)(N downto 0) <= (carry(2)(N+1) & lut_o(3)(N+1) & ccout(2)(N+1 downto 3));
--din(4)(N downto 0) <= (carry(3)(N+1) & lut_o(4)(N+1) & ccout(3)(N+1 downto 3));
--din(5)(N downto 0) <= (carry(4)(N+1) & lut_o(5)(N+1) & ccout(4)(N+1 downto 3));
--din(6)(N downto 0) <= (carry(5)(N+1) & lut_o(6)(N+1) & ccout(5)(N+1 downto 3));
--din(7)(N downto 0) <= (carry(6)(N+1) & lut_o(7)(N+1) & ccout(6)(N+1 downto 3));
--din(8)(N downto 0) <= (carry(7)(N+1) & lut_o(8)(N+1) & ccout(7)(N+1 downto 3));
--din(9)(N downto 0) <= (carry(8)(N+1) & lut_o(9)(N+1) & ccout(8)(N+1 downto 3));
--din(10)(N downto 0) <= (carry(9)(N+1) & lut_o(10)(N+1) & ccout(9)(N+1 downto 3));
--din(11)(N downto 0) <= (carry(10)(N+1) & lut_o(11)(N+1) & ccout(10)(N+1 downto 3));
--din(12)(N downto 0) <= (carry(11)(N+1) & lut_o(12)(N+1) & ccout(11)(N+1 downto 3));
--din(13)(N downto 0) <= (carry(12)(N+1) & lut_o(13)(N+1) & ccout(12)(N+1 downto 3));
--din(14)(N downto 0) <= (carry(13)(N+1) & lut_o(14)(N+1) & ccout(13)(N+1 downto 3));
--din(15)(N downto 0) <= (carry(14)(N+1) & lut_o(15)(N+1) & ccout(14)(N+1 downto 3));

--p <= (ccout(15)(N+1 downto 0) & ccout(14)(1 downto 0) & ccout(13)(1 downto 0) & ccout(12)(1 downto 0) & 
--    ccout(11)(1 downto 0) & ccout(10)(1 downto 0) & ccout(9)(1 downto 0) & ccout(8)(1 downto 0) & 
--    ccout(7)(1 downto 0) & ccout(6)(1 downto 0) & ccout(5)(1 downto 0) & ccout(4)(1 downto 0) & 
--    ccout(3)(1 downto 0) & ccout(2)(1 downto 0) & ccout(1)(1 downto 0) & ccout(0)(1 downto 0));
------------------------------------------------------------------------------------------

end Behavioral;