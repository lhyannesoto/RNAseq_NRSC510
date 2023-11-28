library('gapminder')
gapminder
library('dplyr')

#make a data fram containing columns year, lifeExp, country

select(gapminder, year, lifeExp, country)

select(gapminder, continent, everything())

#order by year using arrange() function
arrange(gapminder, year)

#in descending order
arrange(gapminder, desc(year))

#order by year and lifeExp
arrange (gapminder, year, lifeExp)

#piping cmd+shift+m
#nesting function calls can be hard to read
arrange(select(gapminder, year, lifeExp, country), year, lifeExp)

#now with piping!!!
gapminder %>%
  select(year, lifeExp, country) %>%
  arrange(year, lifeExp)

#FILTERING - taking a subset of my data

#population of >100 000 000
gapminder %>%
  filter(pop>100000000)
#only from asia
gapminder %>%
  filter(pop>100000000) %>%
  filter(continent=='Asia')
#from Brazil or China
gapminder %>%
  filter(pop>100000000) %>%
  filter(country=="Brazil" | country == "China")

##mutate
#make a new column named 'GDP' that equals to multiplying GDP per capita with population
gapminder %>%
  mutate("GDP" = gdpPercap*pop)

#make a new column named 'GDP_bill' that is GDP in billions
gapminder %>%
  mutate("GDP_bill" = gdpPercap * pop/1E9)

#make a new column called cc that pastes the country name followed by the continent, separated by a comma (hint: use the paste function with the sep= argument)
gapminder %>%
  mutate("cc"=paste(country, continent, sep=","))

##SUMMARIZE
#useful for stats: creates a new column but condenses the whole dataset to a number
#let's compute the mean and standard dev of life exp
gapminder %>%
  group_by(continent, year) %>%
  summarize(mu=mean(lifeExp)),
sigma = sd(lifeExp)

#????

#GROUP BY
#group by continent and year
gapminder %>%
  group_by(continent, year)

#4. what's min life exp for each continent and each year


gapminder %>%
  group_by(continent, year) %>%
  summarize(min_life = min(lifeExp)) %>%
  arrange(min_life) #arrange by min life exp

#calculate growth in pop since the first year on record for each country
#a convenient function for you: 'dplyr::first()'

gapminder %>%
  group_by(country) %>%
  arrange(year) %>%
  mutate(rel_growth = pop-first(pop))




















