----------------------------------------------------------------------------------
-- Company: 
-- Engineer: 
-- 
-- Create Date: 08/21/2021 07:51:22 AM
-- Design Name: 
-- Module Name: approximate_multiplier - Behavioral
-- Project Name: 
-- Target Devices: 
-- Tool Versions: 
-- Description: 
-- 
-- Dependencies: 
-- 
-- Revision:
-- Revision 0.01 - File Created
-- Additional Comments:
-- 
----------------------------------------------------------------------------------


library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

-- Uncomment the following library declaration if using
-- arithmetic functions with Signed or Unsigned values
--use IEEE.NUMERIC_STD.ALL;

-- Uncomment the following library declaration if instantiating
-- any Xilinx leaf cells in this code.
library UNISIM;
use UNISIM.VComponents.all;

entity approximate_multiplier is
generic (N : integer:=4); --N:width, M: heigth
port(
    a : in std_logic_vector(N-1 downto 0); --width of the multiplier
    b : in std_logic_vector(N-1 downto 0); --heigth of the multiplier
    d : in std_logic_vector(N downto 0); --has always to be zero, only necessary since generic generate
    p : out std_logic_vector(2*N-1 downto 0)
);
end approximate_multiplier;

architecture Behavioral of approximate_multiplier is

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

--LUT types

constant LUT_A_CON: bit_vector :=X"0000007000000088";
constant LUT_B_CON: bit_vector :=X"0000007000000088";
constant LUT_C_CON: bit_vector :=X"0000007000000088";




--Init signals here

constant PP_NUM_ROWS: integer:= N;
constant PP_NUM_COLS: integer:= N;

type TROW_SIG is array (N downto 0) of std_logic_vector(N+6 downto 0);
signal t_row:TROW_SIG;

type LUT_OUT_SIG is array (N-1 downto 0) of std_logic_vector(N+5 downto 0);
signal lut_out_0:LUT_OUT_SIG;
signal lut_out_1:LUT_OUT_SIG;

type RES_SIG is array (N-1 downto 0) of std_logic_vector(N+5 downto 0);
signal result:RES_SIG;

type CARRY_SIG is array (N-1 downto 0) of std_logic_vector(N+5 downto 0);
signal carry:CARRY_SIG;

--type PP_INT is array (PP_ROWS downto 0) of std_logic_vector(N downto 0);
--signal lut_out: PP_INT; 

begin
--Generate LUTS
PP_ROW_GENERATE: for R in 0 to PP_NUM_ROWS-1 generate --spans multiplier rows
    FIRST_ROW: if R < PP_NUM_ROWS-1 generate
        PP_COL_GENERATE_F1: for C1 in 0 to PP_NUM_COLS-1 generate  --spans multiplier cols
            LUT_A:if C1 < PP_NUM_COLS-1 generate 
                           LUTA00X : LUT6_2  generic map(
                                INIT => LUT_A_CON
                            )
                            port map(
                                I0 => lut_out_1 (R)(C1),
                                I1 => '0',
                                I2 => '0',
                                I3 => t_row(R)(C1),
                                I4 => a(C1),--multiplicand bit
                                I5 => b(R),  --multiplier bit
                                O5 => lut_out_1 (R+1)(C1), --PP_in output
                                O6 => lut_out_0 (R)(C1) --S out put
                            );
                 end generate LUT_A;
--            LUT_B_1:if C1 = PP_NUM_COLS-1 and R=1 generate
--                            LUTB01X : LUT6_2  generic map(
--                                INIT => LUT_B_CON
--                            )
--                            port map(
--                                I0 => '0',
--                                I1 => '0',
--                                I2 => '1',
--                                I3 => t_row(R)(C1),
--                                I4 => a(C1),--multiplicand bit
--                                I5 => b(R),  --multiplier bit
--                                O5 => lut_out_1 (R)(C1), --PP_in output
--                                O6 => lut_out_0 (R)(C1) --S out put
--                            );
--                  end generate LUT_B_1;
--            LUT_B_0:if C1 = PP_NUM_COLS-1 and (R/=1) generate
--                  --LUT B (0)
--                            LUTB00X : LUT6_2  generic map(
--                            INIT => LUT_B_CON
--                            )
--                            port map(
--                                I0 => '0',
--                                I1 => '0',
--                                I2 => '0',
--                                I3 => t_row(R)(C1),
--                                I4 => a(C1),--multiplicand bit
--                                I5 => b(R),  --multiplier bit
--                                O5 => lut_out_1 (R)(C1), --PP_in output
--                                O6 => lut_out_0 (R)(C1) --S out put
--                            );
--                  end generate LUT_B_0;
        end generate PP_COL_GENERATE_F1; --cols1         
     end generate FIRST_ROW;   
--    LAST_ROW: if R = PP_NUM_ROWS-1 generate
--        PP_COL_GENERATE_F2: for C2 in 0 to PP_NUM_COLS generate  --spans multiplier cols
--            LUT_B_0_L:if C2 < PP_NUM_COLS-1 generate 
--                  --LUT_B here
--                     LUTB0LX : LUT6_2  generic map(
--                            INIT => LUT_B_CON
--                            )
--                            port map(
--                                I0 => '0',
--                                I1 => '0',
--                                I2 => '0',
--                                I3 => t_row(R)(C2),
--                                I4 => a(C2),--multiplicand bit
--                                I5 => b(R),  --multiplier bit
--                                O5 => lut_out_1 (R)(C2), --PP_in output
--                                O6 => lut_out_0 (R)(C2) --S out put
--                            );
--                 end generate LUT_B_0_L;
--            LUT_A_L:if C2 = PP_NUM_COLS-1 generate
--                  --LUT A
--                        LUTA0LX : LUT6_2  generic map(
--                                INIT => LUT_A_CON
--                                )
--                                port map(
--                                    I0 => '0',
--                                    I1 => '0',
--                                    I2 => '0',
--                                    I3 => t_row(R)(C2),
--                                    I4 => a(C2),--multiplicand bit
--                                    I5 => b(R),  --multiplier bit
--                                    O5 => lut_out_1 (R)(C2), --PP_in output
--                                    O6 => lut_out_0 (R)(C2) --S out put
--                                );
--                         end generate LUT_A_L;
--            LUT_C_L:if C2 = PP_NUM_COLS generate
--                  --LUT C
--                  LUTC0LX : LUT6_2  generic map(
--                                INIT => LUT_C_CON
--                                )
--                                port map(
--                                    I0 => '0',
--                                    I1 => '0',
--                                    I2 => '0',
--                                    I3 => '0',
--                                    I4 => '0',--multiplicand bit
--                                    I5 => '0',  --multiplier bit
--                                    O5 => lut_out_1 (R)(C2), --PP_in output
--                                    O6 => lut_out_0 (R)(C2) --S out put
--                                );
--                   end generate LUT_C_L;
--        end generate PP_COL_GENERATE_F2; --cols2
--     end generate LAST_ROW;   
   end generate PP_ROW_GENERATE; --rows            
        
--        LUTA00X : LUT6_2
--            generic map(
--                INIT => X"0000007000000088"
--            )
--            port map(
--                I0 => a(0),
--                I1 => '0',
--                I2 => '0',
--                I3 => '0',
--                I4 => '0',
--                I5 => '1'
--              --  O5 => cout(stage)(nr),
--               -- O6 => carry(stage)(nr)
--            );
        
--Generate CARRY4



--PP_ROW_GENERATE_C: for RC in 0 to PP_NUM_ROWS-1 generate --spans multiplier rows
--    FIRST_ROW_CC:if RC < PP_NUM_COLS-1 generate 
--                CC_1: for CC in 0 to (PP_NUM_COLS/4-1) generate --spans multiplier rows
--                            CARRYYX : CARRY4
--                            port map(
--                                CO => carry(RC)(CC*4+4 downto CC*4+1), --only co(3) signal is required
--                                --O => result(RC)(CC*4+3 downto CC*4),
--                                O => t_row(RC)(CC*4+4 downto CC*4+1),
--                                CI => carry(RC)(CC*4), --attach this to co(3) of previous signal
--                                CYINIT => '0',
--                                DI => lut_out_1(RC)(CC*4+3 downto CC*4) ,
--                                S => lut_out_0(RC)(CC*4+3 downto CC*4)                
--                            );
--                      end generate CC_1;
--    end generate FIRST_ROW_CC;
--    LAST_ROW_CC:if RC=PP_NUM_ROWS-1 generate
--                CC_2: for CC1 in 0 to (PP_NUM_COLS/4) generate --spans multiplier rows
--                            CARRYYX_1 : CARRY4
--                            port map(
--                                CO => carry(RC)(CC1*4+4 downto CC1*4+1), --only co(3) signal is required
--                                O => t_row(RC+1)(CC1*4+4 downto CC1*4+1),
--                                CI => carry(RC)(CC1*4), --attach this to co(3) of previous signal
--                                CYINIT => '0',
--                                DI => lut_out_1(RC)(CC1*4+3 downto CC1*4) ,
--                                S => lut_out_0(RC)(CC1*4+3 downto CC1*4)                
--                            );
             
--                end generate CC_2;
--    end generate LAST_ROW_CC;
--end generate PP_ROW_GENERATE_C ;
--    CARRYYX : CARRY4
--    port map(
--        CO => carry(RC) (nr*4+4 downto nr*4+1), --only co(3) signal is required
--        O => result(RC)(),
--        CI => carry(RC) (nr*4), --attach this to co(3) of previous signal
--        CYINIT => '0',
--        DI => lut_out_1(RC)(CC1 downto CC1) ,
--        S => lut_out_0(RC)()                
--    );

--p<=(t_row(PP_NUM_ROWS)(PP_NUM_COLS downto 0) & "0000000");


end Behavioral;