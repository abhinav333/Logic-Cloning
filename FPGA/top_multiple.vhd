library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

entity top_multiple is
generic (N : integer :=16; COP: integer:=1);
port(
topclk : in std_logic;
topam : in std_logic_vector(N-1 downto 0);
topbm : in std_logic_vector(N-1 downto 0);
topdm : in std_logic_vector(N downto 0);
toppm : out std_logic_vector(2*N*COP-1 downto 0)
);
end top_multiple;

architecture Behavioral of top_multiple is


component top 
generic (N : integer ; M : integer );
port(
topclk : in std_logic;
topa : in std_logic_vector(N-1 downto 0);
topb : in std_logic_vector(M-1 downto 0);
topd : in std_logic_vector(N downto 0);
topp : out std_logic_vector(N+M-1 downto 0)
);
end component;

begin

gen_cop: for i in 0 to COP-1 generate
mul_wrapper : top
generic map(
    N => N,
    M => N
)
port map(
    topclk => topclk,
    topa => topam,
    topb => topbm,
    topd => topdm,
    topp => toppm((i+1)*2*N-1 downto i*2*N)
);
end generate;
end Behavioral;
