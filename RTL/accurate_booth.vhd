----------------------------------------------------------------------------------
-- Company: 
-- Engineer: Abhinav
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
generic (N : integer); --N:width, M: heigth
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

constant LUT_A1_CON: bit_vector :=X"C66CC66C0DDF0880";
constant LUT_A_CON: bit_vector :=X"F0C3C3A55A3C3CF0";
constant LUT_A2_CON: bit_vector :=X"3666999333333333";
constant LUT_B_CON: bit_vector :=X"0F3C3C3CC3C3C30F";

type TROW_SIG is array (N/2 downto 0) of std_logic_vector(N+5 downto 0);
signal t_row:TROW_SIG:= (others => (others => '0'));

type TROW_SIG1 is array (N/2 downto 0) of std_logic_vector(1 downto 0);
signal t_row_1:TROW_SIG1:= (others => (others => '0'));

type DI_SIG is array (N/2-1 downto 0) of std_logic_vector(N+5 downto 0);
signal d_row:DI_SIG:= (others => (others => '0'));


type CARRY_SIG is array (N/2-1 downto 0) of std_logic_vector(N+5 downto 0);
signal carry:CARRY_SIG:= (others => (others => '0'));


type LUT_SIG is array (N/2-1 downto 0) of std_logic_vector(N+5 downto 0);
signal lut_out_0:LUT_SIG:= (others => (others => '0'));

type LUT_SIG1 is array (N/2-1 downto 0) of std_logic_vector(1 downto 0);
signal lut_out_1:LUT_SIG1:= (others => (others => '0'));


signal pp_prod:std_logic_vector(N-1 downto 0):=  (others => '0');

signal b_app:std_logic_vector(N downto 0):=  (others => '0');



begin
--Generate LUTS
t_row_1 (0)(1 downto 0)<="11";
b_app<=b & '0';
PP_ROW_GENERATE: for R in 0 to (N/2-1) generate --spans multiplier rows
    PP_COL_GENERATE_F1: for C1 in 0 to (N+1) generate  --spans multiplier cols
                FIRST_LUT: if C1=0 generate
                    LUTA00X : LUT6_2  generic map(
                                    INIT => LUT_A1_CON
                                )
                                port map(
                                    I0 => a(C1),
                                    I1 => t_row(R)(C1+1),
                                    I2 => b_app(2*R),
                                    I3 => b_app(2*R+1),--multiplicand bit
                                    I4 => b_app(2*R+2),  --multiplier bit
                                    I5 => '1',
                                    O5 => lut_out_1(R)(0), --PP_in output
                                    O6 => lut_out_0(R)(C1) --S output
                                );
                    pp_prod(2*R)<=lut_out_0(R)(C1);
                    carry(R)(0)<=lut_out_1(R)(0);        
                 end generate FIRST_LUT;
                NEXT_LUT:if C1 > 0 and C1 < N generate 
                        LUTA00X : LUT6_2  generic map(
                                INIT => LUT_A_CON
                            )
                            port map(
                                I0 => a(C1-1),
                                I1 => a(C1),
                                I2 => t_row(R)(C1+1),
                                I3 => b_app(2*R),
                                I4 => b_app(2*R+1),--multiplicand bit
                                I5 => b_app(2*R+2),  --multiplier bit
                                O6 => lut_out_0(R)(C1) --S output
                            );
                    d_row(R)(C1-1)<=t_row(R)(C1+1);
                 end generate NEXT_LUT;
                 NEXT_LUT1:if C1 = N generate 
                        LUTA00X : LUT6_2  generic map(
                                INIT => LUT_A2_CON
                            )
                            port map(
                                I0 => a(N-1),
                                I1 => t_row_1(R)(0),
                                I2 => b_app(2*R),
                                I3 => b_app(2*R+1),--multiplicand bit
                                I4 => b_app(2*R+2),  --multiplier bit
                                I5 => '1',
                                O5 => lut_out_1(R)(1), --PP_inoutput
                                O6 => lut_out_0(R)(C1) --S output
                            );
                    d_row(R)(C1-1)<=lut_out_1(R)(1);                   
                 end generate NEXT_LUT1;
         
           NEXT_LUT2:if C1 = N+1 generate 
                        LUTA00X : LUT6_2  generic map(
                                INIT => LUT_B_CON
                            )
                            port map(
                                I0 => a(N-1),
                                I1 => a(N-1),
                                I2 => t_row_1(R)(1),
                                I3 => b_app(2*R),
                                I4 => b_app(2*R+1),--multiplicand bit
                                I5 => b_app(2*R+2),  --multiplier bit
                                O6 => lut_out_0(R)(C1) --S output
                            );
                 d_row(R)(C1-1)<=t_row_1(R)(1);
                 end generate NEXT_LUT2;
                         
           end generate PP_COL_GENERATE_F1;
end generate PP_ROW_GENERATE ;
        
RR_1:for R in 0 to (N/2-1) generate
--d_row(R)(N downto 0)<=t_row_1(R)(1) & lut_out_1(R)(1) & t_row(R)(N downto 2);
CC_1: for CC in 0 to (N/4) generate --spans multiplier rows
                            CARRYYX : CARRY4
                            port map(
                                CO => carry(R)(CC*4+4 downto CC*4+1), --only co(3) signal is required
                                O => t_row(R+1)(CC*4+3 downto CC*4),
                                CI => carry(R)(CC*4), --attach this to co(3) of previous signal
                                CYINIT => '0',
                                DI =>d_row(R)(CC*4+3 downto CC*4),
                                S => lut_out_0(R)(CC*4+4 downto CC*4+1)
                            );
                      end generate CC_1;
                      t_row_1(R+1)(0)<=carry(R)(N+1);
                      t_row_1(R+1)(1)<=carry(R)(N+1);
                      pp_prod(R*2+1)<=t_row(R+1)(0);
    end generate RR_1;
p<=t_row(N/2)(N downto 1) & pp_prod(N-1 downto 0);
end Behavioral;