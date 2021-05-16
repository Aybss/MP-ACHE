#include <stdio.h>

int floatAdd(int firstOper, int secondOper);
unsigned int typespecify(int operand, unsigned int sign, 
                          unsigned int exponent, unsigned int significand);
int noInput(int firstOper, int secondOper, unsigned int firstType, 
                  unsigned int secondType, unsigned int signPart, 
                   unsigned int signPart2);
int unlimitedInp(int firstOper, int secondOper, unsigned int signPart, unsigned int signPart2,
                    unsigned int firstType, unsigned int secondType);
int zeroInputTypes(int firstOper, int secondOper, unsigned int firstType, 
                unsigned int secondType);
int normalData(int firstOper, int secondOper, unsigned int signPart, unsigned int signPart2, 
                    unsigned int firstExpo, unsigned int secondExpo, 
                      unsigned int firstSig, unsigned int secondSig,
                          unsigned int unNormalized);
int lookOVFlow(unsigned int exponent);
int countOvFlow(unsigned int exponent);
int stickPart(unsigned int significand, unsigned int parted);
int floatAdd(int firstOper, int secondOper)
{
  int partOutput;
  unsigned int firstType, secondType, unNormalized;
  unsigned int signPart, firstExpo, firstSig;
  unsigned int signPart2, secondExpo, secondSig;
  #define DO_DEBUGGING_PRINTS 0
  signPart = (firstOper & 0x80000000) >> 31;
  firstExpo = (firstOper & 0x7F800000) >> 23;
  firstSig = ((firstOper << 9) >> 9) & 0x007FFFFF;
  firstType = typespecify(firstOper, signPart, firstExpo, firstSig);
  signPart2 = (secondOper & 0x80000000) >> 31;
  secondExpo = (secondOper & 0x7F800000) >> 23;
  secondSig = ((secondOper << 9) >> 9) & 0x007FFFFF;
  secondType = typespecify(secondOper, signPart2, secondExpo, secondSig);
  unNormalized = 0;
  if(firstType == 2 || secondType == 2)
  {
    partOutput = zeroInputTypes(firstOper, secondOper, firstType, secondType);
  }
  else if(firstType == 1 || secondType == 1)
  {
    partOutput = unlimitedInp(firstOper, secondOper, signPart, signPart2, firstType, secondType);
  } 
  else if(firstType == 0 || secondType == 0)
  {
    partOutput = noInput(firstOper, secondOper, firstType, secondType, signPart, signPart2);
  }
  else if(firstType == 4 && secondType == 4)
  {
    unNormalized = 1;
    firstExpo = firstExpo + 1;
    secondExpo = secondExpo + 1;
    //add
    partOutput = normalData(firstOper, secondOper, signPart, signPart2, firstExpo,
                            secondExpo, firstSig, secondSig, unNormalized);
  }
  if(firstType == 3 && secondType == 4)
  {
    secondExpo = secondExpo + 1;
    firstSig = firstSig | 0x00800000;
    partOutput = normalData(firstOper, secondOper, signPart, signPart2, firstExpo,
                            secondExpo, firstSig, secondSig, unNormalized);
  }
  else if(firstType == 4 && secondType == 3)
  {
    firstExpo = firstExpo + 1;
    secondSig = secondSig | 0x00800000;
    partOutput = normalData(firstOper, secondOper, signPart, signPart2, firstExpo,
                            secondExpo, firstSig, secondSig, unNormalized);
  }
  else if(firstType == 3 && secondType == 3)
  {
    firstSig = firstSig | 0x00800000;
    secondSig = secondSig | 0x00800000;
    partOutput = normalData(firstOper, secondOper, signPart, signPart2, firstExpo,
                            secondExpo, firstSig, secondSig, unNormalized);
  }
  #if DO_DEBUGGING_PRINTS
    fprintf(stderr, "-----Operand 1-----\n");
    fprintf(stderr, "Sign: %01x\n", signPart);
    fprintf(stderr, "Exponent: %08x\n", firstExpo);
    fprintf(stderr, "Significand: %08x\n", firstSig);
  #endif
  #if DO_DEBUGGING_PRINTS
    fprintf(stderr, "-----Operand 2-----\n");
    fprintf(stderr, "Sign: %01x\n", signPart2);
    fprintf(stderr, "Exponent: %08x\n", secondExpo);
    fprintf(stderr, "Significand: %08x\n", secondSig);
  #endif
 
  return partOutput;
}
unsigned int typespecify(int operand, unsigned int sign, 
                          unsigned int exponent, unsigned int significand)
{
  unsigned int secondDel = 5;
  if(exponent == 0x00)
  {
    if(significand == 0x000000)
    {
      secondDel = 0;
      if(sign == 0)
      {
        #if DO_DEBUGGING_PRINTS
          fprintf(stderr, "Value %08x = +0\n", operand);
        #endif
      }
      else if(sign == 1)
      {
        #if DO_DEBUGGING_PRINTS
          fprintf(stderr, "Value %08x = -0\n", operand);
        #endif
      }
    }
    else
    {
      secondDel = 4;
      #if DO_DEBUGGING_PRINTS
        fprintf(stderr, "Value %08x is De-Normalized\n", operand);
      #endif
    }
  }
  else if(exponent == 0xFF)
  {
    if(significand == 0x000000)
    {
      secondDel = 1;
      if(sign == 0)
      {
        #if DO_DEBUGGING_PRINTS
          fprintf(stderr, "Value %08x = +Infinity\n", operand);
        #endif
      }
      else if(sign == 1)
      {
        #if DO_DEBUGGING_PRINTS
          fprintf(stderr, "Value %08x = -Infinity\n", operand);
        #endif
      }
    }
    else
    {
      secondDel = 2;
      #if DO_DEBUGGING_PRINTS
        fprintf(stderr, "Value %08x is NaN\n", operand);
      #endif
    }
  }
  else
  {
      secondDel = 3;
      #if DO_DEBUGGING_PRINTS
        fprintf(stderr, "Value %08x is Normalized\n", operand);
      #endif
  }
  return secondDel;
}
int noInput(int firstOper, int secondOper, unsigned int firstType, 
                unsigned int secondType, unsigned int signPart, 
                  unsigned int signPart2)
{
  int secondDel;
  unsigned int sign;
  if(firstType == 0 && secondType != 0)
  {
    secondDel = secondOper;
  }
  else if(secondType == 0 && firstType != 0)
  {
    secondDel = firstOper;
  }
  else if(secondType == 0 && firstType == 0)
  {
    if(signPart < signPart2)
    {
      sign = signPart;
    }
    else if(signPart > signPart2)
    {
      sign = signPart2;
    }
    else
    {
      sign = signPart;
    }
    secondDel = ((firstOper & 0x7FFFFFFF) | (sign << 31));
  }

  return secondDel;
}
int unlimitedInp(int firstOper, int secondOper, unsigned int signPart, unsigned int signPart2,
                    unsigned int firstType, unsigned int secondType)
{
  int secondDel;
  if(firstType == 1 && secondType != 1)
  {
    secondDel = firstOper;
  }
  else if(secondType == 1 && firstType != 1)
  {
    secondDel = secondOper;
  }
  else if(secondType == 1 && firstType == 1)
  {
    if(signPart == signPart2)
    {
      secondDel = firstOper;
    }
    else
    {
      secondDel = 0xFFC00000;
    }
  }
  return secondDel;
}
int zeroInputTypes(int firstOper, int secondOper, unsigned int firstType, 
                  unsigned int secondType)
{
  int secondDel;
  if(firstType == 2 && secondType != 2)
  {
    secondDel = (firstOper | 0x7FC00000);
  }
  else if(secondType == 2 && firstType != 2)
  {
    secondDel = (secondOper | 0x7FC00000);
  }
  else if(secondType == 2 && firstType == 2)
  {
    secondDel = (secondOper | 0x7FC00000);
  }
  return secondDel;
}
int normalData(int firstOper, int secondOper, unsigned int signPart, unsigned int signPart2, 
                    unsigned int firstExpo, unsigned int secondExpo, 
                      unsigned int firstSig, unsigned int secondSig,
                          unsigned int unNormalized)
{
  unsigned int sign, exponent, significand, ogsignificand; 
  unsigned int parted, firstDel, secondDel, pointDecimal, minused;
  unsigned int caseMinus = 0;
  if(firstExpo > secondExpo)
    {
      sign = signPart;
      firstDel = firstExpo - secondExpo;
      significand = ((secondSig << 2) & 0x03FFFFFC);
      ogsignificand = ((firstSig << 2) & 0x03FFFFFC);
      exponent = secondExpo;
    }
    else if(firstExpo < secondExpo)
    {
      sign = signPart2;
      firstDel = secondExpo - firstExpo;
      significand = ((firstSig << 2) & 0x03FFFFFC);
      ogsignificand = ((secondSig << 2) & 0x03FFFFFC);
      exponent = firstExpo;
    }
    else
    {
      firstDel = 0;
      if(firstSig < secondSig)
      {
        sign = signPart2;
        significand = ((firstSig << 2) & 0x03FFFFFC);
        ogsignificand = ((secondSig << 2) & 0x03FFFFFC);
      }
      else if(firstSig > secondSig)
      {
        sign = signPart;
        significand = ((secondSig << 2) & 0x03FFFFFC);
        ogsignificand = ((firstSig << 2) & 0x03FFFFFC);
      }
      else
      {
        if(signPart == signPart2)
        {
          sign = signPart;
          significand = ((secondSig << 2) & 0x03FFFFFC);
          ogsignificand = ((firstSig << 2) & 0x03FFFFFC);
        }
        else
        {
          sign = 0;
          significand = 0x00000000;
          ogsignificand = 0x00000000;
          caseMinus = 1;
        }
      }
      exponent = firstExpo;
    }
    if(firstDel < 0)
    {
      firstDel = -firstDel;
    }
    parted = 0;
    pointDecimal = 0;
    int i = 0;
    while(i < firstDel)
    {
      if((significand & 0x00000001) == 0x00000001)
      {
        parted = 1;
      }
      significand = significand >> 1;
      exponent++;
      pointDecimal = stickPart(significand, parted);
      i++;
    }
    if(signPart == signPart2)
    {
      significand = ogsignificand + significand;
      minused = 0;      
    }
    else
    {
      significand = ogsignificand - significand;
      minused = 1;
    }
    pointDecimal = stickPart(significand, parted);
    if(significand == 0x00000000)
    {
      secondDel = 0x00000000;
    }
    else
    {      
      if((significand & 0xFFFFFFFF) < 0x02000000) 
      {
        while((significand & 0xFFFFFFFF) < 0x02000000)
        {
          if(unNormalized == 1 && exponent <= 0x01)
          {
            exponent = 0x00;
            caseMinus = 1;
            break;
          }
          significand = significand << 1;
          exponent--;
          pointDecimal = stickPart(significand, parted);
          if (lookOVFlow(exponent) == 1)
          {
            exponent = 0x00;
            significand = 0x00000000;
            caseMinus = 1;
            break;
          }
        }
      }
      if((significand & 0xFFFFFFFF) > 0x03FFFFFF) 
      {
        while((significand & 0xFFFFFFFF) > 0x03FFFFFF)
        {
          if((significand & 0x00000001) == 0x00000001)
          {
            parted = 1;
          }
          significand = ((significand >> 1) & 0x7FFFFFFF);
          exponent++;
          pointDecimal = stickPart(significand, parted);
          if (countOvFlow(exponent) == 1)
          {
            exponent = 0xFF;
            significand = 0x00000000;
            caseMinus = 1;
            break;
          }
        }
      }
      significand = ((significand >> 2) & 0x00FFFFFF);
      if(caseMinus != 1)
      { 
        if(pointDecimal == 2)
        {
          if((significand & 0x00000001) == 0x00000001)
          {
            significand = significand + 1;
          }
        }
        else if((pointDecimal == 3) && (minused == 0))
        {
          significand = significand + 1;
        }
        else if(pointDecimal == 4)
        {
          significand = significand + 1;
        }
        if((significand & 0xFFFFFFFF) > 0x00FFFFFF) 
        {
          while((significand & 0xFFFFFFFF) > 0x00FFFFFF)
          {
            significand = ((significand >> 1) & 0x7FFFFFFF);
            exponent++;
            if (countOvFlow(exponent) == 1)
            {
              exponent = 0xFF;
              significand = 0x00000000;
              caseMinus = 1;
              break;
            }
          }
        }
      }
      significand = (significand & 0x007FFFFF);
      secondDel = (significand | ((exponent << 23) & 0x7F800000));
      secondDel = ((secondDel & 0x7FFFFFFF) | ((sign & 0x00000001) << 31));
    }
  
  return secondDel;
}
int lookOVFlow(unsigned int exponent)
{
  unsigned int secondDel = 0;
  if(exponent > 0x800000FE)
  {
    secondDel = 1;
  }
  return secondDel;
}
int countOvFlow(unsigned int exponent)
{
  unsigned int secondDel = 0;
  if(exponent >= 0x000000FF)
  {
    secondDel = 1;
  }
  return secondDel;
}
int stickPart(unsigned int significand, unsigned int parted)
{
  unsigned int floatPoint, secondDel;
  floatPoint = (significand & 0x00000003);
  secondDel = floatPoint + parted;
  if((secondDel == 2) && (floatPoint == 1))
  {
    secondDel = 1;
  }
  else if((secondDel == 3) && (floatPoint == 3))
  {
    secondDel = 4;
  }

  return secondDel;
}
